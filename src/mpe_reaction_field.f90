!****h* FHI-aims/mpe_reaction_field
!  NAME
!    mpe_reaction_field - calculates the reaction field of the MPE implicit
!                         solvent model
!  SYNOPSIS

module mpe_reaction_field

!  PURPOSE
!    This module provides subroutines to evaluate the reaction field.
!    For more information regarding the MPE implicit solvent model
!    please also see:
!     [1] M. Sinstein, C. Scheurer, S. Matera, V. Blum, K. Reuter, 
!           and H. Oberhofer, submitted
!     [2] D. Rinaldi, A. Bouchy, J.L. Rivail, and V. Dillet, 
!           J. Chem. Phys. 120, 2343 (2004)
!     [3] V. Dillet, D. Rinaldi, J.G. Angyan, and J.L. Rivail, 
!           Chem. Phys. Lett. 202, 18 (1993)
!     [4] J.L. Rivail and D. Rinaldi, Chem. Phys. 18, 233-242 (1976)
!     [5] J.G. Kirkwood, J. Chem. Phys. 2, 351-361 (1934)
!  USES

   use mpe_constants, only: MPE_CONST
   use mpe_types, only: &
         Basis, &
         SolHarmBasis, &
         BasisCenter, &
         SolHarmBasisCenter, &
         SolHarmBasisCenter_size_lm, &
         SolHarmBasisCenter_lbound_lm, &
         DielectricContinuum, &
         DielectricInterface, &
         InterfacePoint, &
         InterfacePoint_vector_extract, &
         SpVecTuple, &
         SpVecTuple_vector_from_dense_vector, &
         SpVecTuple_vector_mpi_bcast

   ! The following modules are part of FHI-aims and are thus to be
   ! replaced by corresponding functionalities when the MPE model
   ! is to be completely decoupled from FHI-aims.
   use types, only: dp
   use constants, only: sqrt_pi
   use mpi_tasks, only: f_stop => aims_stop
   use aims_memory_tracking, only: &
         f_allocate => aims_allocate, &
         f_deallocate => aims_deallocate
   use localorb_io, only: localorb_info, OL_norm, use_unit
   use timing, only: get_timestamps, get_times, output_times
   use synchronize_mpi_basic, only: sync_vector
   use cartesian_ylm, only: &
         n_max_cartesian, &
         initialize_cartesian_ylm, &
         evaluate_onecenter_cartesians, &
         tab_ylm_onecenter_cartesian, &
         evaluate_onecenter_cartesian_gradient_terms
   use runtime_choices, only: ifp, &
         charge

   implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   private  ! make everything private by default

   ! PUBLIC ROUTINES (subroutines and functions)

   public   calculate_reaction_field_coefficients
   public   calculate_potential_at_points
   public   initialize_mpe_solver
   public   adjust_basis_functions
   public   cleanup_mpe_solver

   type, private :: MPEFactorizationMemory
      real(dp), allocatable, private :: QR(:,:)
      real(dp), allocatable, private :: QR_tau(:)
      real(dp), allocatable, private :: U(:,:)
      real(dp), allocatable, private :: S(:)
      real(dp), allocatable, private :: VT(:,:)
   end type

   ! The Settings type encapsulates all information that needs to be
   ! stored during subsequent calls to the MPE solver routines.
   ! It is meant to be instanced by the caller which has to use the
   ! initialization (initialize_mpe_solver) and finalization
   ! (cleanup_mpe_solver) routines.
   type, public :: MPESettings
      logical, public :: use_scalapack = .false.
      integer, public :: factorization_type = MPE_CONST % FACTZN_UNDEF
      real(dp), public :: sparsity_threshold = 0.e0_dp

      ! BLACS context
      integer, private :: mpe_communicator
      integer, private :: mpe_blacs_ctxt
      integer, private :: mpe_setup_blacs_ctxt

      ! storage for iterative SLE solver,
      ! only allocated if iterative solver is activated
      type(MPEFactorizationMemory), private :: factorizations

      ! certain initialization only has to be done once, thus remember
      logical, private ::  module_initialized = .false.
      logical, private ::  solver_first_run = .true.
   end type


   abstract interface
      subroutine CBHartreePot(n, points, v_h, v_h_gradient)
         use types, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: points(3,n)
         real(dp), intent(out) :: v_h(n)
         real(dp), intent(out), optional :: v_h_gradient(3,n)
      end subroutine
   end interface

   ! EXTERNAL FUNCTIONS

   integer, external :: PILAENV, NUMROC, INDXG2L, INDXG2P, INDXL2G, ICEIL
   real(dp), external :: DNRM2, DLAMCH, PDLAMCH


   ! PRIVATE PARAMETERS

   ! constant terms
   real(dp), parameter :: ONE = 1.e0_dp, ZERO = 0.e0_dp

   ! ScaLAPACK descriptor
   integer, parameter :: DLEN_ = 9
   integer, parameter :: BLOCK_CYCLIC_2D = 1
   integer, parameter :: DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, &
                         NB_ = 6, RSRC_ = 7, CSRC_ = 8, LLD_ = 9

   integer :: sc_desc_U(DLEN_), sc_desc_VT(DLEN_)

contains

!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/calculate_reaction_field_coefficients
!  NAME
!    calculate_reaction_field_coefficients
!  SYNOPSIS

subroutine calculate_reaction_field_coefficients( settings, &
      geometry_changed, continua, interfaces, v_hartree_callback, adjR2 )

!  PURPOSE
!    This is just a wrapper routine around the actual calculation routines.
!
!  USES
   implicit none
!  ARGUMENTS
   integer, parameter :: RHS_cols = 1

   type(MPESettings), intent(inout) :: settings
   logical, intent(in) :: geometry_changed
   type(DielectricContinuum), intent(inout) :: continua(:)
   type(DielectricInterface), intent(in) :: interfaces(:)
   real(dp), intent(out) :: adjR2(RHS_cols)
   type(DielectricInterface), allocatable :: iftemp(:)
   logical, allocatable :: mask(:)
   integer :: i_if, i_iftemp
   procedure(CBHartreePot) :: v_hartree_callback

   intrinsic count
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o geometry_changed -- did the problem's geometry change since last call?
!   o continua -- dielectric continua with assigned basis functions
!   o interfaces -- dielectric interfaces
!   o v_hartree_callback -- callback routine evaluating the Hartree potential
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o continua -- dielectric continua, now including expansion coefficients
!   o adjR2 -- coefficient of determination of the solved SLE
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'calculate_reaction_field_coefficients'
   character(132) :: info_str

   if (.not. settings % module_initialized) then
      write(info_str,'(A)') 'Module has not been initialized.'
      call localorb_info(info_str)
      call f_stop('Module has not been initialized', func)
   endif

   allocate(mask(size(interfaces)))
   do i_if = 1, size(interfaces), 1
      mask(i_if) = interfaces(i_if)%n_bc .ne. 0
   end do
   allocate(iftemp(count(mask)))
   i_iftemp = 1
   do i_if = 1, size(interfaces), 1
      if (mask(i_if)) then
         iftemp(i_iftemp) = interfaces(i_if)
         i_iftemp = i_iftemp + 1
      endif
   end do
   deallocate(mask)

   if (settings % use_scalapack) then
      call calculate_reaction_field_coefficients_scalapack ( &
               settings = settings, &
               geometry_changed = geometry_changed, &
               continua = continua, &
               interfaces = iftemp, &
               v_hartree_callback = v_hartree_callback, &
               adjR2 = adjR2 )
   else
      call calculate_reaction_field_coefficients_lapack ( &
               settings = settings, &
               geometry_changed = geometry_changed, &
               continua = continua, &
               interfaces = iftemp, &
               v_hartree_callback = v_hartree_callback, &
               adjR2 = adjR2 )
   endif

   deallocate(iftemp)

   settings % solver_first_run = .false.

end subroutine calculate_reaction_field_coefficients
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/initialize_mpe_reaction_field
!  NAME
!    initialize_mpe_reaction_field
!  SYNOPSIS

subroutine initialize_mpe_solver( settings, mpi_comm )

!  PURPOSE
!    This subroutine checks whether MPE parameters have been specified in the
!    control.in and assigns default values if this is not the case.
!
!  USES
   implicit none
!  ARGUMENTS
   type(MPESettings), intent(inout) :: settings
   integer, intent(in) :: mpi_comm
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o mpi_comm -- MPI communication handle to be used by this module
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'initialize_mpe_reaction_field'

   integer :: myid, np, npcol, nprow

   if (settings % module_initialized) &
      call cleanup_mpe_solver( settings )

   settings % mpe_communicator = mpi_comm

   ! BLACS GRID INIT
   if (settings % use_scalapack) then
      ! get number of processors
      call BLACS_PINFO(myid, np)
      ! use all processes in a grid closest to quadratic shape
      do npcol = nint(sqrt(real(np,dp))), 2, -1
         if (mod(np,npcol).eq.0 ) exit
      enddo
      ! at the end of the above loop, np is always divisible by npcol
      nprow = np/npcol
      ! initialize BLACS grids
      settings % mpe_blacs_ctxt = mpi_comm
      call BLACS_GRIDINIT( settings % mpe_blacs_ctxt, &
                           'Column-major', nprow, npcol )
      settings % mpe_setup_blacs_ctxt = mpi_comm
      call BLACS_GRIDINIT( settings % mpe_setup_blacs_ctxt, &
                           'Column-major', np, 1 )
   endif ! use_scalapack

   settings % module_initialized = .true.
   settings % solver_first_run = .true.

end subroutine initialize_mpe_solver
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/cleanup_mpe_solver
!  NAME
!    cleanup_mpe_solver
!  SYNOPSIS

subroutine cleanup_mpe_solver( settings )

!  PURPOSE
!    This subroutine deallocates all persistent MPE variables except
!    the reaction field on the integration grid which is deallocated
!    in another routine.
!
!  USES
   implicit none
!  ARGUMENTS
   type(MPESettings), intent(inout) :: settings
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   if (.not. settings % module_initialized) return

   if (allocated(settings % factorizations % QR)) &
      call f_deallocate(settings % factorizations % QR, &
            name="+QR factorized Left-Hand Side")
   if (allocated(settings % factorizations % QR_tau)) &
      call f_deallocate(settings % factorizations % QR_tau, &
            name="Scalar factors for QR factorization")
   if (allocated(settings % factorizations % U)) &
      call f_deallocate(settings % factorizations % U, &
            name="+SVD: Left Singular Vectors")
   if (allocated(settings % factorizations % S)) &
      call f_deallocate(settings % factorizations % S, &
            name="SVD: Singular Values")
   if (allocated(settings % factorizations % VT)) &
      call f_deallocate(settings % factorizations % VT, &
            name="+SVD: Right Singular Vectors")

   ! BLACS GRID EXIT
   if (settings % use_scalapack) then
      call BLACS_GRIDEXIT(settings % mpe_blacs_ctxt)
      call BLACS_GRIDEXIT(settings % mpe_setup_blacs_ctxt)
   endif ! use_scalapack

   settings % module_initialized = .false.
   settings % solver_first_run = .true.

end subroutine cleanup_mpe_solver
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/adjust_basis_functions
!  NAME
!    adjust_basis_functions
!  SYNOPSIS

subroutine adjust_basis_functions(continua, interfaces)

!  PURPOSE
!    This subroutine allows to adjust the basis functions to a new
!    geometry of dielectric interfaces.
!
!  USES
   implicit none
!  ARGUMENTS
   type(DielectricContinuum), intent(inout) :: continua(:)
   type(DielectricInterface), intent(in) :: interfaces(:)
!  INPUTS
!   o continua -- dielectric continua with assigned basis functions
!   o interfaces -- dielectric interfaces
!  OUTPUT
!   o continua -- dielectric continua with adjusted basis functions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   integer :: i_dc, i_if, i_center, i_point
   real(dp) :: relvec(3), min_distsq, max_distsq, distsq, new_rscale
   character(132) :: info_str

   do i_dc = lbound(continua,1), ubound(continua,1)
      select type (b => continua(i_dc)%basis)
      class is (SolHarmBasis)

         do i_center = lbound(b%centers,1), ubound(b%centers,1)

            min_distsq = huge(ONE)
            max_distsq = tiny(ONE)
            ! only consider points in contact with the continuum to which the basis belongs
            do i_if = 1, size(interfaces)
               if ( (i_dc.eq.interfaces(i_if)%dc_ind_pos .or. &
                     i_dc.eq.interfaces(i_if)%dc_ind_neg) .and. &
                    interfaces(i_if)%n_bc .ne. 0 ) then

                  do i_point = lbound(interfaces(i_if)%p,1), &
                                 ubound(interfaces(i_if)%p,1)
                     relvec = interfaces(i_if)%p(i_point)%coord - &
                              b%centers(i_center)%coord
                     distsq = dot_product(relvec, relvec)
                     min_distsq = min(min_distsq, distsq)
                     max_distsq = max(max_distsq, distsq)
                  enddo
               endif
            enddo ! i_if

            ! pick correct scaling factor
            new_rscale = ONE
            if ( min_distsq.ne.huge(ONE) .and. max_distsq.ne.tiny(ONE) ) then
               select case(b%solharm_type)
                  case(MPE_CONST%BASIS_REG)
                     new_rscale = ONE/sqrt(max_distsq)
                  case(MPE_CONST%BASIS_IRR)
                     new_rscale = ONE/sqrt(min_distsq)
               endselect
            endif

            ! re-scaling of coefficients is not necessary as long as the 
            !  coefficients are formulated in terms of the unsaled coordinates

            ! set new scaling
            b%centers(i_center)%rscale = new_rscale

         enddo ! i_center

      class default
         call f_stop('Internal Error: adjustment for this basis is '//&
                        'not implemented.')
      end select
   enddo ! i_dc

end subroutine adjust_basis_functions
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/calculate_reaction_field_coefficients_lapack
!  NAME
!    calculate_reaction_field_coefficients_lapack
!  SYNOPSIS

subroutine calculate_reaction_field_coefficients_lapack( &
      settings, geometry_changed, continua, interfaces, v_hartree_callback, adjR2 )

!  PURPOSE
!    Set up and solve a system of linear equations in order to
!    calculate the reaction field coefficients R_l'm'
!
!  USES
   implicit none

!  ARGUMENTS
   integer, parameter :: RHS_cols = 1

   type(MPESettings), intent(inout) :: settings
   logical, intent(in) :: geometry_changed
   type(DielectricContinuum), intent(inout) :: continua(:)
   type(DielectricInterface), intent(in) :: interfaces(:)
   real(dp), intent(out) :: adjR2(RHS_cols)   ! adjusted R^2
   procedure(CBHartreePot) :: v_hartree_callback
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o geometry_changed -- did the geometry change since the last call ?
!                         (i.e. is a different LHS matrix is to be expected?)
!   o continua -- dielectric continua with assigned basis functions
!   o interfaces -- dielectric interfaces
!   o v_hartree_callback -- callback routine evaluating the Hartree potential
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o continua -- contains now also the expansion coefficients of the basis
!   o adjR2 -- coefficient of determination of the solved SLE
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

!HOOKLAPACK
   character(*), parameter :: func = 'calculate_reaction_field_coefficients'//&
                                       '_lapack'
   character(132) :: info_str
   integer :: info

   ! variables for timings
   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm
   real(dp) :: cpu_time, clock_time

   ! switches
   logical :: do_setup_LHS
   logical :: perform_factorization
   logical :: need_qr
   logical :: need_svd
   logical :: svd_on_qr = .true.
   logical :: iterative_solver_converged

   ! left-hand side matrix LHS and right-hand side matrix RHS in the system of
   ! linear equations LHS * x = RHS
   ! x is the matrix of reaction field coefficients
   real(dp), allocatable :: LHS(:,:), RHS(:,:)

   ! dimensions and sizes
   integer, allocatable :: n_rows_if(:), n_cols_dc(:), n_rows_dc(:)
   integer :: n_variables
   integer :: i_if, n_interfaces, i_dc, n_continua
   integer :: LHS_rows, LHS_cols, RHS_rows
   integer, parameter :: IND_RF = 1 !TODO replace

   ! variables for statistical analysis
   real(dp) :: TSS(RHS_cols)     ! variation of Y
   real(dp) :: RSS(RHS_cols)     ! sum of squared residuals
   real(dp) :: RMSD(RHS_cols)    ! root mean squared deviation
   real(dp) :: R2(RHS_cols)      ! coefficient of determination

   integer :: myid
   integer :: i_center, i_sv

   ! get process id
   call MPI_Comm_rank(settings % mpe_communicator, myid, info)

   ! DIMENSIONS
   write(info_str,'(4X,A)') 'determine problem dimensions'
   call localorb_info(info_str)

   ! define dimensions of SLE
   n_interfaces = size(interfaces)
   n_continua = size(continua)
   allocate(n_rows_if(n_interfaces))
   allocate(n_cols_dc(n_continua))

   ! The columns consist of one block for each continuum basis
   do i_dc = 1, n_continua, 1
      n_cols_dc(i_dc) = continua(i_dc)%basis%get_size()
   enddo ! i_dc

   ! The rows consist of one block for each interface
   do i_if = 1, n_interfaces, 1
      n_rows_if(i_if) = interfaces(i_if)%n_bc * size(interfaces(i_if)%p)
   enddo ! i_if

   ! In interface case, charge conservation needs to be enforced
   if (ifp) then
      allocate(n_rows_dc(n_continua))
      do i_dc = 1, n_continua, 1
         select type(b => continua(i_dc)%basis)
            class is (SolHarmBasis)
               select case(b%solharm_type)
                  case(MPE_CONST%BASIS_REG)
                     n_rows_dc(i_dc) = 0
                  case(MPE_CONST%BASIS_IRR)
                     ! Currently, this is only total charge
                     n_rows_dc(i_dc) = 1
                  case default
                     call f_stop('Internal Error: MPE with interface plane and '//&
                             'given basis type not implemented')
               end select ! b%solharm_type
            class default
               call f_stop('Internal Error: MPE with interface plane and '//&
                       'given basis type not implemented')
         end select ! type(b => continua(i_dc)%basis)
      enddo ! i_dc
   endif ! ifp

   ! grade of determination
   write(info_str,'(6X,A2,2(2X,A12),2X,A8)') 'interface', '#conditions', &
                                             '#variables', 'ratio'
   call localorb_info(info_str)
   do i_if = 1, n_interfaces
      n_variables = n_cols_dc(interfaces(i_if)%dc_ind_pos) + &
                        n_cols_dc(interfaces(i_if)%dc_ind_neg)
      write(info_str,'(6X,I2,2(2X,I12),2X,F8.3)') i_if, n_rows_if(i_if), &
                        n_variables, n_rows_if(i_if)/real(n_variables,dp)
      call localorb_info(info_str)
   enddo

   LHS_rows = sum(n_rows_if)
   if (ifp) &
           LHS_rows = LHS_rows + sum(n_rows_dc)
   LHS_cols = sum(n_cols_dc)
   RHS_rows = max(LHS_rows, LHS_cols) ! TODO: This doesn't seem right

   write(info_str,'(6X,A,I10)') 'total number of conditions: ', LHS_rows
   call localorb_info(info_str)
   write(info_str,'(6X,A,I10)') 'total number of variables:  ', LHS_cols
   call localorb_info(info_str)
   write(info_str,'(8X,A,F8.3)') 'ratio conds/vars:  ', &
                           real(LHS_rows,dp) / real(LHS_cols,dp)
   call localorb_info(info_str)

   ! reset this
   adjR2(:) = 0.e0_dp
   TSS(:)   = 0.e0_dp
   RSS(:)   = 0.e0_dp
   RMSD(:)  = 0.e0_dp
   R2(:)    = 0.e0_dp

   ! MAIN PART
   if (myid.eq.0) then

      ! SET UP MATRICES
      do_setup_LHS = settings % solver_first_run .or. geometry_changed
      if (do_setup_LHS) then
         write(info_str,'(4X,A)') 'setting up matrices'
      else
         write(info_str,'(4X,A)') 'setting up right-hand side'
      endif
      call localorb_info(info_str)

      ! get timestamp
      call get_timestamps(cpu_time, clock_time)

      if (do_setup_LHS) then
         ! calculate left-hand side matrix
         call f_allocate(LHS, LHS_rows, LHS_cols, name="+Left-Hand Side")
         LHS = 0.e0_dp
         do i_if = 1, n_interfaces
            call serial_setup_LHS_part ( interface = interfaces(i_if), &
                     continuum = continua(interfaces(i_if)%dc_ind_pos), &
                     positive_side = .true., &
                     LHS = LHS, &
                     first_global_row = 1+sum(n_rows_if(1:i_if-1)), &
                     first_global_col = &
                        1+sum(n_cols_dc(1:interfaces(i_if)%dc_ind_pos-1)) )
            call serial_setup_LHS_part ( interface = interfaces(i_if), &
                     continuum = continua(interfaces(i_if)%dc_ind_neg), &
                     positive_side = .false., &
                     LHS = LHS, &
                     first_global_row = 1+sum(n_rows_if(1:i_if-1)), &
                     first_global_col = &
                        1+sum(n_cols_dc(1:interfaces(i_if)%dc_ind_neg-1)) )
         enddo ! i_if

         ! Internal conditions at the end
         if (ifp) then
            do i_dc = 1, n_continua, 1
               call serial_setup_internal_LHS_part(continuum = continua(i_dc), &
                        LHS = LHS, &
                        first_global_row = &
                                    1+sum(n_rows_if)+sum(n_rows_dc(1:i_dc-1)), &
                        first_global_col = &
                                    1+sum(n_cols_dc(1:i_dc-1)) )
            enddo ! i_dc
         endif ! ifp
      endif ! do_setup_LHS

      ! calculate right-hand side vector
      call f_allocate(RHS, RHS_rows, RHS_cols, name="+Right-Hand Side")
      RHS = 0.e0_dp
      do i_if = 1, n_interfaces
         call serial_setup_RHS_part( &
               interface = interfaces(i_if), &
               delta_inv_eps = &
                  1.e0_dp/continua(interfaces(i_if)%dc_ind_neg)%eps - &
                  1.e0_dp/continua(interfaces(i_if)%dc_ind_pos)%eps, &
               v_hartree_callback = v_hartree_callback, &
               RHS = RHS, &
               first_global_row = 1+sum(n_rows_if(1:i_if-1)) )
      enddo ! i_if

      ! Internal conditions at the end
      !
      ! TODO: This is a bit inconsistent, as serial_setup_internal_LHS_part
      !       works only on one continuum, but serial_setup_internal_RHS_part
      !       on all at once
      if (ifp) then
         call serial_setup_internal_RHS_part(continua = continua, &
                 RHS = RHS, &
                 weights = (/ 0.e0_dp, .5e0_dp, .5e0_dp /), &
                 first_global_row = 1+sum(n_rows_if) )
         ! TODO: figure out how to get explicit multipole components in there,
         !       if conditions for higher multipoles ever get implemented
      endif ! ifp

      ! sum up time
      call get_times(cpu_time, clock_time, unsynced=.true.)
      call output_times(deffmt, 'Setting up matrices', &
                      cpu_time, clock_time, defprio)

      ! STATISTICAL ANALYSIS 1
      call serial_statistical_analysis_1( RHS, TSS )

      !MS: This is a dummy, iterative solver goes here
      iterative_solver_converged = .false.
      perform_factorization = do_setup_LHS .and. &
            (.not. iterative_solver_converged)
      need_qr = (settings % factorization_type.eq.MPE_CONST % FACTZN_QR) &
         .or.(settings % factorization_type.eq.MPE_CONST % FACTZN_QRpSVD)
      need_svd = (settings % factorization_type.eq.MPE_CONST % FACTZN_SVD) &
         .or.(settings % factorization_type.eq.MPE_CONST % FACTZN_QRpSVD)

      if ( perform_factorization ) then
         write(info_str,'(4X,A)') 'factorizing coefficient matrix'
         call localorb_info(info_str)

         if ( need_qr ) then
            call get_timestamps(cpu_time, clock_time)
            ! TODO: This is just a hotfix, should actually check if number
            ! of boundary points has changed
            if (geometry_changed .and. &
                    allocated(settings % factorizations % QR)) then
               call f_deallocate(settings % factorizations % QR, &
                     name="+QR factorized Left-Hand Side")
            endif

            if (.not. allocated(settings % factorizations % QR) ) then
               call f_allocate(settings % factorizations % QR, &
                     LHS_rows, LHS_cols, name="+QR factorized Left-Hand Side")
            endif
            if (.not. allocated(settings % factorizations % QR_tau) ) then
               call f_allocate(settings % factorizations % QR_tau, &
                     LHS_cols, name="Scalar factors for QR factorization")
            endif

            call serial_factorization_QR( &
                  LHS, &
                  settings % factorizations % QR, &
                  settings % factorizations % QR_tau )

            call get_times(cpu_time, clock_time, unsynced=.true.)
            call output_times(deffmt, 'QR factorization', cpu_time, clock_time, defprio)

         else
            svd_on_qr = .false.

         endif ! need_qr

         if ( need_svd ) then
            call get_timestamps(cpu_time, clock_time)
            if (.not. allocated(settings % factorizations % S) ) then
               call f_allocate(settings % factorizations % S, &
                     LHS_cols, name="SVD: Singular Values")
            endif
            if (.not. allocated(settings % factorizations % VT) ) then
               call f_allocate(settings % factorizations % VT, &
                     LHS_cols, LHS_cols, name="+SVD: Right Singular Vectors")
            endif
            if ( svd_on_qr ) then
               if (.not. allocated(settings % factorizations % U) ) then
                  call f_allocate(settings % factorizations % U, &
                     LHS_cols, LHS_cols, name="+SVD: Left Singular Vectors")
               endif
               ! copy R from QR decomposition into LHS
               call DLASET("Lower", LHS_cols-1, LHS_cols-1, &
                     ZERO, ZERO, &
                     LHS(2,1), LHS_rows )
               call DLACPY("Upper", LHS_cols, LHS_cols, &
                     settings % factorizations % QR, LHS_rows, &
                     LHS, LHS_rows )

               call serial_factorization_SVD( &
                     LHS, &
                     settings % factorizations % U, &
                     settings % factorizations % S, &
                     settings % factorizations % VT, &
                     use_driver = 'custom  ')
            else
               ! If SVD is preformed on full LHS, shape might have changed
               ! TODO: This is just a hotfix, should actually check if number
               ! of boundary points has changed
               if (geometry_changed .and. &
                       allocated(settings % factorizations % U)) then
                  call f_deallocate(settings % factorizations % U, &
                        name="+SVD: Left Singular Vectors")
               endif

               if (.not. allocated(settings % factorizations % U) ) then
                  call f_allocate(settings % factorizations % U, &
                        LHS_rows, LHS_cols, name="+SVD: Left Singular Vectors")
               endif

               call serial_factorization_SVD( &
                     LHS, &
                     settings % factorizations % U, &
                     settings % factorizations % S, &
                     settings % factorizations % VT, &
                     use_driver = 'SVD     ' )
            endif ! have qr_factorization

            call get_times(cpu_time, clock_time, unsynced=.true.)
            call output_times(deffmt, 'SVD factorization', cpu_time, clock_time, defprio)
         endif ! need_svd
      endif ! perform_factorization

      if (allocated(LHS)) &
         call f_deallocate(LHS, name="+Left-Hand Side")

      if (.not. iterative_solver_converged) then
         write(info_str,'(4X,A)') 'direct solution'
         call localorb_info(info_str)
         ! get timestamp
         call get_timestamps(cpu_time, clock_time)

         if ( need_qr ) then
            call serial_solver_Q( &
                  settings % factorizations % QR, &
                  settings % factorizations % QR_tau, &
                  RHS, RSS )
         else
            RSS = ZERO
         endif
         if ( need_svd ) then
            call serial_solver_SVD( &
                  settings % factorizations % U, &
                  settings % factorizations % S, &
                  settings % factorizations % VT, &
                  RHS, RSS )
         else
            call serial_solver_R( &
                  settings % factorizations % QR, &
                  RHS, RSS ) !TODO: actually do something with the RSS here
         endif

         call get_times(cpu_time, clock_time, unsynced=.true.)
         call output_times(deffmt, 'Direct solution', cpu_time, clock_time, defprio)
      endif ! .not. iterative_solver_converged

      ! STATISTICAL ANALYSIS 2
      call statistical_analysis_2( LHS_rows, LHS_cols, &
            TSS, RSS, RMSD, R2, adjR2 )

      ! EXTRACTION
      ! extract reaction field coefficients from RHS matrix,
      ! which has been overwritten by the solution matrix
      call serial_extract_rfc( RHS, continua, settings % sparsity_threshold )

      if (allocated(RHS)) &
         call f_deallocate(RHS, name="+Right-Hand Side")

   endif ! myid.eq.0


   ! SYNCHRONIZATION

   ! get timestamp
   call get_timestamps(cpu_time, clock_time)

   !TODO move to subroutine
   ! sync only expansion coefficients,
   ! all the statistics is entirely done on first thread

   ! adjR2 was only calculated on 0th proc
   call sync_vector(adjR2, RHS_cols)

   do i_dc = 1, n_continua
      select type (b => continua(i_dc)%basis)
      class is (SolHarmBasis)
         do i_center = lbound(b%centers,1), &
                       ubound(b%centers,1)
            call SpVecTuple_vector_mpi_bcast( b%centers(i_center)%coeff, &
                        0, settings % mpe_communicator )
         enddo ! i_center
      class default
         call f_stop('Internal Error: synchronization for this basis is '//&
                        'not implemented.')
      end select
   enddo ! i_dc

   ! sum up time
   call get_times(cpu_time, clock_time)
   call output_times(deffmt, 'Synchronization', &
             cpu_time, clock_time, defprio)


   ! ANALYSIS

   ! get timestamp
   call get_timestamps(cpu_time, clock_time)

   ! STATISTICAL ANALYSIS
   call print_statistical_analysis_rfc( use_unit, &
                           TSS(1), RSS(1), RMSD(1), R2(1), adjR2(1) )

   ! SPARSITY ANALYSIS
   call print_sparsity_analysis_rfc( use_unit, continua(IND_RF)%basis )

   ! sum up time
   call get_times(cpu_time, clock_time)
   call output_times(deffmt, 'Time for fit analysis', &
             cpu_time, clock_time, defprio)

end subroutine calculate_reaction_field_coefficients_lapack
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/calculate_reaction_field_coefficients_scalapack
!  NAME
!    calculate_reaction_field_coefficients_scalapack
!  SYNOPSIS

subroutine calculate_reaction_field_coefficients_scalapack( &
      settings, geometry_changed, continua, interfaces, v_hartree_callback, adjR2 )

!  PURPOSE
!    Set up and solve a system of linear equations in order to
!    calculate the reaction field coefficients R_l'm'
!
!  USES
   implicit none

!  ARGUMENTS
   integer, parameter :: glo_RHS_cols = 1

   type(MPESettings), intent(inout) :: settings
   logical, intent(in) :: geometry_changed
   type(DielectricContinuum), intent(inout) :: continua(:)
   type(DielectricInterface), intent(in) :: interfaces(:)
   real(dp), intent(out) :: adjR2(glo_RHS_cols)   ! adjusted R^2
   procedure(CBHartreePot) :: v_hartree_callback
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o geometry_changed -- did the geometry change since the last call ?
!                         (i.e. is a different LHS matrix is to be expected?)
!   o continua -- dielectric continua with assigned basis functions
!   o interfaces -- dielectric interfaces
!   o v_hartree_callback -- callback routine evaluating the Hartree potential
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!   o continua -- contains now also the expansion coefficients of the basis
!   o adjR2 -- coefficient of determination of the solved SLE
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

!HOOKSCALAPACK
   character(*), parameter :: func = 'calculate_reaction_field_coefficients'//&
                                       '_scalapack'
   character(132) :: info_str
   integer :: info

   ! variables for timings
   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm
   real(dp) :: cpu_time, clock_time

   ! switches
   logical :: do_setup_LHS
   logical :: perform_factorization
   logical :: need_qr
   logical :: need_svd
   logical :: svd_on_qr = .true.
   logical :: iterative_solver_converged

   ! left-hand side matrix LHS and right-hand side matrix RHS in the system of
   ! linear equations LHS * x = RHS
   real(dp), allocatable :: loc_LHS(:,:), loc_RHS(:,:)
   integer :: loc_LHS_rows, loc_LHS_cols, sc_desc_LHS(DLEN_), &
              loc_RHS_rows, loc_RHS_cols, sc_desc_RHS(DLEN_)

   ! temporary variables used for LHS/RHS matrix set-up
   real(dp), allocatable :: loc_setup_LHS(:,:), loc_setup_RHS(:,:)
   integer :: loc_setup_LHS_rows, loc_setup_LHS_cols, sc_desc_setup_LHS(DLEN_), &
              loc_setup_RHS_rows, loc_setup_RHS_cols, sc_desc_setup_RHS(DLEN_)

   ! temporary variables used for LHS/RHS matrix set-up, internal conditions
   real(dp), allocatable :: loc_setup_LHS_int(:,:), loc_setup_RHS_int(:,:)
   integer :: loc_setup_LHS_rows_int, loc_setup_LHS_cols_int, sc_desc_setup_LHS_int(DLEN_), &
              loc_setup_RHS_rows_int, loc_setup_RHS_cols_int, sc_desc_setup_RHS_int(DLEN_), &
              glo_setup_LHS_rows_int, glo_setup_RHS_rows_int, &
              glo_setup_LHS_rows_int_padded, glo_setup_RHS_rows_int_padded

   ! reduced system of LHS/RHS matrices
   integer :: loc_red_LHS_rows, loc_red_LHS_cols, sc_desc_red_LHS(DLEN_), &
              loc_red_RHS_rows, loc_red_RHS_cols, sc_desc_red_RHS(DLEN_)

   ! descriptors for SVD
   integer :: sc_desc_red_LHS_in_full(DLEN_)

   ! dimensions and sizes
   integer, allocatable :: glo_n_rows_if(:), glo_n_cols_dc(:)
   integer, allocatable :: glo_n_rows_dc(:)
   integer :: n_variables
   integer :: i_if, n_interfaces, i_dc, n_continua
   integer :: glo_LHS_rows,     glo_LHS_cols,     glo_RHS_rows
   integer :: glo_LHS_rows_tot,                   glo_RHS_rows_tot
   integer, parameter :: IND_RF = 1 !TODO replace

   ! variables for statistical analysis
   real(dp) :: TSS(glo_RHS_cols)     ! variation of Y
   real(dp) :: RSS(glo_RHS_cols)     ! sum of squared residuals
   real(dp) :: RMSD(glo_RHS_cols)    ! root mean squared deviation
   real(dp) :: R2(glo_RHS_cols)      ! coefficient of determination

   integer :: myid
   integer :: block_divisor
   integer :: BLACS_MYPNUM, BLACS_NPROCS, irow

   ! DEBUG: coordinate scaling with interface
   integer :: i_sv

   ! get process id
   call MPI_Comm_rank(settings % mpe_communicator, myid, info)

   ! DIMENSIONS
   write(info_str,'(4X,A)') 'determine problem dimensions'
   call localorb_info(info_str)

   ! define dimensions of SLE
   n_interfaces = size(interfaces)
   n_continua = size(continua)
   allocate(glo_n_rows_if(n_interfaces))
   allocate(glo_n_cols_dc(n_continua))

   ! The columns consist of one block for each continuum basis
   do i_dc = 1, n_continua, 1
      glo_n_cols_dc(i_dc) = continua(i_dc)%basis%get_size()
   enddo ! i_dc

   ! The rows consist of one block for each interface
   do i_if = 1, n_interfaces, 1
      glo_n_rows_if(i_if) = interfaces(i_if)%n_bc * size(interfaces(i_if)%p)
   enddo ! i_if

   ! In interface case, charge conservation needs to be enforced
   if (ifp) then
      allocate(glo_n_rows_dc(n_continua))
      do i_dc = 1, n_continua, 1
         select type(b => continua(i_dc)%basis)
            class is (SolHarmBasis)
               select case(b%solharm_type)
                  case(MPE_CONST%BASIS_REG)
                     glo_n_rows_dc(i_dc) = 0
                  case(MPE_CONST%BASIS_IRR)
                     ! Currently, this is only total charge
                     glo_n_rows_dc(i_dc) = 1
                  case default
                     call f_stop('Internal Error: MPE with interface plane and '//&
                                     'given basis type not implemented')
               end select ! b%solharm_type
            class default
               call f_stop('Internal Error: MPE with interface plane and '//&
                               'given basis type not implemented')
         end select ! type(b => continua(i_dc)%basis)
      enddo ! i_dc
   endif ! ifp

   ! grade of determination
   write(info_str,'(6X,A2,2(2X,A12),2X,A8)') 'interface', '#conditions', &
                                             '#variables', 'ratio'
   call localorb_info(info_str)
   do i_if = 1, n_interfaces
      n_variables = glo_n_cols_dc(interfaces(i_if)%dc_ind_pos) + &
                        glo_n_cols_dc(interfaces(i_if)%dc_ind_neg)
      write(info_str,'(6X,I2,2(2X,I12),2X,F8.3)') i_if, glo_n_rows_if(i_if), &
                        n_variables, glo_n_rows_if(i_if)/real(n_variables,dp)
      call localorb_info(info_str)
   enddo

   glo_LHS_rows = sum(glo_n_rows_if)
   if (ifp) then
           glo_setup_LHS_rows_int = sum(glo_n_rows_dc)
           glo_setup_RHS_rows_int = glo_setup_LHS_rows_int

           ! TODO:
           ! This is a dirty, dirty trick. Basically, this sets up an internal
           ! condition matrix way larger than necessary, but ultimately only the
           ! first two rows are filled and copied. Reason: DESCINIT in
           ! get_scalapack_descriptors_for_SLE_setup would otherwise give a warning
           ! because loc_*_rows would be 0 on all processes except the 0th.
           call BLACS_PINFO(BLACS_MYPNUM, BLACS_NPROCS)
           glo_setup_LHS_rows_int_padded = glo_setup_LHS_rows_int * BLACS_NPROCS
           glo_setup_RHS_rows_int_padded = glo_setup_RHS_rows_int * BLACS_NPROCS


           glo_LHS_rows_tot = glo_LHS_rows + glo_setup_LHS_rows_int
   else
           glo_LHS_rows_tot = glo_LHS_rows
   endif

   glo_LHS_cols = sum(glo_n_cols_dc)

   glo_RHS_rows     = max(glo_LHS_rows    , glo_LHS_cols)
   glo_RHS_rows_tot = max(glo_LHS_rows_tot, glo_LHS_cols) ! TODO: look into this again

   write(info_str,'(6X,A,I10)') 'total number of conditions: ', glo_LHS_rows_tot
   call localorb_info(info_str)
   write(info_str,'(6X,A,I10)') 'total number of variables:  ', glo_LHS_cols
   call localorb_info(info_str)
   write(info_str,'(6X,A,F8.3)') 'ratio conds/vars:  ', &
                           real(glo_LHS_rows_tot,dp) / real(glo_LHS_cols,dp)
   call localorb_info(info_str)

   ! reset this
   adjR2(:) = 0.e0_dp
   TSS(:)   = 0.e0_dp
   RSS(:)   = 0.e0_dp
   RMSD(:)  = 0.e0_dp
   R2(:)    = 0.e0_dp

   ! get descriptors for full system, setup
   block_divisor = product(interfaces(:)%n_bc)
   call get_scalapack_descriptors_for_SLE_setup( &
                  settings % mpe_setup_blacs_ctxt, &
                  glo_LHS_rows, glo_LHS_cols, &
                  glo_RHS_rows, glo_RHS_cols, &
                  block_divisor, &
                  loc_setup_LHS_rows, loc_setup_LHS_cols, sc_desc_setup_LHS, &
                  loc_setup_RHS_rows, loc_setup_RHS_cols, sc_desc_setup_RHS )

   if (ifp) then
      call get_scalapack_descriptors_for_SLE_setup( &
                  settings % mpe_setup_blacs_ctxt, &
                  glo_setup_LHS_rows_int_padded, glo_LHS_cols, &
                  glo_setup_RHS_rows_int_padded, glo_RHS_cols, &
                  glo_setup_LHS_rows_int, & ! TODO: currently moves all relevant elements to 0th process
                  loc_setup_LHS_rows_int, loc_setup_LHS_cols_int, sc_desc_setup_LHS_int, &
                  loc_setup_RHS_rows_int, loc_setup_RHS_cols_int, sc_desc_setup_RHS_int)
   endif

   ! get descriptors for full system, calculation
   call get_scalapack_descriptors_for_SLE( &
                  settings % mpe_blacs_ctxt, &
                  glo_LHS_rows_tot, glo_LHS_cols, &
                  glo_RHS_rows_tot, glo_RHS_cols, &
                  loc_LHS_rows, loc_LHS_cols, sc_desc_LHS, &
                  loc_RHS_rows, loc_RHS_cols, sc_desc_RHS )

   ! get descriptors for reduced system, calculation
   call get_scalapack_descriptors_for_SLE( &
                  settings % mpe_blacs_ctxt, &
                  glo_LHS_cols, glo_LHS_cols, &
                  glo_LHS_cols, glo_RHS_cols, &
                  loc_red_LHS_rows, loc_red_LHS_cols, sc_desc_red_LHS, &
                  loc_red_RHS_rows, loc_red_RHS_cols, sc_desc_red_RHS )

   write(info_str,'(6X,A,I4,A,I4)') 'matrix block size:  ', &
                        sc_desc_LHS(NB_), ' x ', sc_desc_LHS(MB_)
   call localorb_info(info_str)


   ! SET UP MATRICES
   do_setup_LHS = settings % solver_first_run .or. geometry_changed
   if (do_setup_LHS) then
      write(info_str,'(4X,A)') 'setting up matrices'
   else
      write(info_str,'(4X,A)') 'setting up right-hand side'
   endif
   call localorb_info(info_str)

   ! get timestamp
   call get_timestamps(cpu_time, clock_time)

   ! calculate left-hand side matrix
   if (do_setup_LHS) then
      call f_allocate(loc_LHS, loc_LHS_rows, loc_LHS_cols, name="+Left-Hand Side")
      call f_allocate(loc_setup_LHS, loc_setup_LHS_rows, loc_setup_LHS_cols, &
            name="+Setup of Left-Hand Side")
      loc_LHS = 0.e0_dp
      loc_setup_LHS = 0.e0_dp
      do i_if = 1, n_interfaces
         call parallel_setup_local_LHS_part ( interface = interfaces(i_if), &
                  continuum = continua(interfaces(i_if)%dc_ind_pos), &
                  positive_side = .true., &
                  part = loc_setup_LHS, sc_desc=sc_desc_setup_LHS, &
                  first_global_row = 1+sum(glo_n_rows_if(1:i_if-1)), &
                  first_global_col = &
                     1+sum(glo_n_cols_dc(1:interfaces(i_if)%dc_ind_pos-1)) )
         call parallel_setup_local_LHS_part ( interface = interfaces(i_if), &
                  continuum = continua(interfaces(i_if)%dc_ind_neg), &
                  positive_side = .false., &
                  part = loc_setup_LHS, sc_desc=sc_desc_setup_LHS, &
                  first_global_row = 1+sum(glo_n_rows_if(1:i_if-1)), &
                  first_global_col = &
                     1+sum(glo_n_cols_dc(1:interfaces(i_if)%dc_ind_neg-1)) )
      enddo ! i_if
      ! re-distribute matrix for better computational performance
      call PDGEMR2D(glo_LHS_rows, glo_LHS_cols, &
                  loc_setup_LHS, 1, 1, sc_desc_setup_LHS, &
                  loc_LHS, 1, 1, sc_desc_LHS, sc_desc_LHS(CTXT_) )
      call f_deallocate(loc_setup_LHS, name="+Setup of Left-Hand Side")

      ! Internal conditions at the end
      if (ifp) then
         if ( (myid .eq. 0 .and. loc_setup_LHS_rows_int .ne. glo_setup_LHS_rows_int) ) &
            call f_stop("Internal error: Internal conditions must all be set up on 0th process "// &
                    "in the current implementation.")

         call f_allocate(loc_setup_LHS_int, loc_setup_LHS_rows_int, loc_setup_LHS_cols_int, &
                 name="+Setup of Left-Hand Side, internal conditions")
         loc_setup_LHS_int = 1.e0_dp ! yes, 1. v.i.

         ! Foolproofing the matrix copying stuff. This should NEVER actually play a role
         ! Here, check if anything has been written yet where the internal conditions
         ! will go.
         call PDGEMR2D(glo_setup_LHS_rows_int, glo_LHS_cols, &
                 loc_LHS, glo_LHS_rows+1, 1, sc_desc_LHS, &
                 loc_setup_LHS_int, 1, 1, sc_desc_setup_LHS_int, sc_desc_setup_LHS_int(CTXT_))
         ! 0th process should now contain the zeros from the redistributed matrix
         ! All other processes should still contain ones from after allocation
         if ( (myid .eq. 0 .and. any(loc_setup_LHS_int .ne. 0.e0_dp)) .or. &
              (myid .ne. 0 .and. any(loc_setup_LHS_int .ne. 1.e0_dp)) ) &
            call f_stop("Internal error: Internal block of LHS matrix already contains elements, "//&
                        "or top rows of setup matrix are not on 0th process")
         ! Re-initialize from 0, just to be sure
         loc_setup_LHS_int = 0.e0_dp

         do i_dc = 1, n_continua, 1
            if (myid .eq. 0) then
               if (loc_setup_LHS_cols_int .ne. glo_LHS_cols) &
                  call f_stop("Internal error: Serial subroutine will not work if matrix is distributed over "//&
                     "more than 1 pcol")
               call serial_setup_internal_LHS_part(continuum = continua(i_dc), &
                        LHS = loc_setup_LHS_int, &
                        first_global_row = 1+sum(glo_n_rows_dc(1:i_dc-1)), &
                        first_global_col = 1+sum(glo_n_cols_dc(1:i_dc-1)) )
            endif ! myid == 0
         enddo ! i_dc
         ! re-distribute matrix for better computational performance
         call PDGEMR2D(glo_setup_LHS_rows_int, glo_LHS_cols, &
                     loc_setup_LHS_int, 1, 1, sc_desc_setup_LHS_int, &
                     loc_LHS, glo_LHS_rows+1, 1, sc_desc_LHS, sc_desc_LHS(CTXT_) )

         ! More foolproofing. Again, should never play a role.
         ! Here, check if the right rows (the top 2, which should always be on the
         ! 0th process) have been copied
         loc_setup_LHS_int = 0.e0_dp
         call PDGEMR2D(glo_setup_LHS_rows_int, glo_LHS_cols, &
                 loc_LHS, glo_LHS_rows+1, 1, sc_desc_LHS, &
                 loc_setup_LHS_int, 1, 1, sc_desc_setup_LHS_int, sc_desc_setup_LHS_int(CTXT_))
         ! On 0th process: Filled rows should have been copied, overwritten by
         ! zeros, and copied back. No row should be all zeros, except it contains
         ! only centers with lmin>0, which should not happen unless someone hacks
         ! the code.
         !
         ! On all other processes: Nothing should have been copied.
         do irow = lbound(loc_setup_LHS_int, dim=1), ubound(loc_setup_LHS_int, dim=1)
            if ( (myid .eq. 0 .and. all(loc_setup_LHS_int(irow,:) .eq. 0.e0_dp)) .or. &
                 (myid .ne. 0 .and. any(loc_setup_LHS_int(irow,:) .ne. 0.e0_dp)) ) then
               write(info_str, '(A)') " "
               call localorb_info(info_str)
               write(info_str, '(A)') " * Internal error: It seems like the wrong rows for internal conditions "//&
                           "have been copied to the MPE matrix."
               call localorb_info(info_str)
               write(info_str, '(A)') "                   This error will also appear if all centers "//&
                           "for an external potential have lmin>0."
               call localorb_info(info_str)
               write(info_str, '(A)') "                   The latter should however only happen if you "//&
                           "modified the code."
               call localorb_info(info_str)
               write(info_str, '(A)') "                   You can turn off this sanity check at your own risk."
               call localorb_info(info_str)
               call f_stop("Internal error: It seems like the wrong rows for internal conditions "//&
                           "have been copied to the MPE matrix.")
            endif
         enddo ! irow

         call f_deallocate(loc_setup_LHS_int, name="+Setup of Left-Hand Side, internal conditions")
      endif ! ifp

   endif ! do_setup_LHS

   ! calculate right-hand side vector
   call f_allocate(loc_RHS, loc_RHS_rows, loc_RHS_cols, name="+Right-Hand Side")
   call f_allocate(loc_setup_RHS, loc_setup_RHS_rows, loc_setup_RHS_cols, &
         name="+Setup of Right-Hand Side")
   loc_RHS = 0.e0_dp
   loc_setup_RHS = 0.e0_dp
   do i_if = 1, n_interfaces
      call parallel_setup_local_RHS_part( &
            interface = interfaces(i_if), &
            delta_inv_eps = &
               1.e0_dp/continua(interfaces(i_if)%dc_ind_neg)%eps - &
               1.e0_dp/continua(interfaces(i_if)%dc_ind_pos)%eps, &
            v_hartree_callback = v_hartree_callback, &
            loc_RHS = loc_setup_RHS, sc_desc=sc_desc_setup_RHS, &
            first_global_row = 1+sum(glo_n_rows_if(1:i_if-1)) )
   enddo ! i_if
   ! re-distribute matrix for better computational performance
   call PDGEMR2D(glo_RHS_rows, glo_RHS_cols, &
               loc_setup_RHS, 1, 1, sc_desc_setup_RHS, &
               loc_RHS, 1, 1, sc_desc_RHS, sc_desc_RHS(CTXT_) )
   call f_deallocate(loc_setup_RHS, name="+Setup of Right-Hand Side")

   ! Internal conditions at the end
   if (ifp) then
      if ( (myid .eq. 0 .and. loc_setup_RHS_rows_int .ne. glo_setup_RHS_rows_int) ) &
         call f_stop("Internal error: Internal conditions must all be set up on 0th process "// &
                 "in the current implementation.")
      if ( loc_setup_RHS_cols_int .ne. 1 ) &
         call f_stop("Internal error: RHS matrix must have exactly 1 column")

      call f_allocate(loc_setup_RHS_int, loc_setup_RHS_rows_int, loc_setup_RHS_cols_int, &
            name="+Setup of Right-Hand Side, internal conditions")
      loc_setup_RHS_int = 1.e0_dp ! yes, 1. v.i.

      ! In analogy to above:
      ! Foolproofing the matrix copying stuff. This should NEVER actually play a role
      ! Here, check if anything has been written yet where the internal conditions
      ! will go.
      call PDGEMR2D(glo_setup_RHS_rows_int, glo_RHS_cols, &
              loc_RHS, glo_RHS_rows+1, 1, sc_desc_RHS, &
              loc_setup_RHS_int, 1, 1, sc_desc_setup_RHS_int, sc_desc_setup_RHS_int(CTXT_))
      ! 0th process should now contain the zeros from the redistributed matrix
      ! All other processes should still contain ones from after allocation
      if ( (myid .eq. 0 .and. any(loc_setup_RHS_int .ne. 0.e0_dp)) .or. &
           (myid .ne. 0 .and. any(loc_setup_RHS_int .ne. 1.e0_dp)) ) &
         call f_stop("Internal error: Internal block of RHS matrix already contains elements, "//&
                     "or top rows of setup matrix are not on 0th process")
      ! Re-initialize from 0, just to be sure
      loc_setup_RHS_int = 0.e0_dp

      ! TODO: This is a bit inconsistent, as serial_setup_internal_LHS_part
      !       works only on one continuum, but serial_setup_internal_RHS_part
      !       on all at once
      if (myid .eq. 0) then
         call serial_setup_internal_RHS_part(continua = continua, &
                 RHS = loc_setup_RHS_int, &
                 weights = (/ 0.e0_dp, .5e0_dp, .5e0_dp /), &
                 first_global_row = 1 )
         ! TODO: figure out how to get explicit multipole components in there,
         !       if conditions for higher multipoles ever get implemented
      endif ! myid == 0

      ! re-distribute matrix for better computational performance
      call PDGEMR2D(glo_setup_RHS_rows_int, glo_RHS_cols, &
                  loc_setup_RHS_int, 1, 1, sc_desc_setup_RHS_int, &
                  loc_RHS, glo_RHS_rows+1, 1, sc_desc_RHS, sc_desc_RHS(CTXT_) )

      ! Unfortunately, no general post-copying sanity check is possible here, because
      ! unlike the LHS, the RHS can (and often will) be zero.

      call f_deallocate(loc_setup_RHS_int, name="+Setup of Right-Hand Side, internal conditions")
   endif ! ifp


   ! sum up time
   call get_times(cpu_time, clock_time)
   call output_times(deffmt, 'Setting up matrices', &
          cpu_time, clock_time, defprio)

   ! STATISTICAL ANALYSIS 1
   call parallel_statistical_analysis_1( loc_RHS, sc_desc_RHS, TSS )

   !MS: This is a dummy, iterative solver goes here
   iterative_solver_converged = .false.
   perform_factorization = do_setup_LHS .and. &
         (.not. iterative_solver_converged)
   need_qr = (settings % factorization_type.eq.MPE_CONST % FACTZN_QR) &
      .or.(settings % factorization_type.eq.MPE_CONST % FACTZN_QRpSVD)
   need_svd = (settings % factorization_type.eq.MPE_CONST % FACTZN_SVD) &
      .or.(settings % factorization_type.eq.MPE_CONST % FACTZN_QRpSVD)

   if ( perform_factorization ) then
      write(info_str,'(4X,A)') 'factorizing coefficient matrix'
      call localorb_info(info_str)

      if ( need_qr ) then
         call get_timestamps(cpu_time, clock_time)
         ! TODO: This is just a hotfix, should actually check if number
         ! of boundary points has changed
         if (geometry_changed .and. &
                 allocated(settings % factorizations % QR)) then
            call f_deallocate(settings % factorizations % QR, &
                  name="+QR factorized Left-Hand Side")
         endif

         if (.not. allocated(settings % factorizations % QR) ) then
            call f_allocate(settings % factorizations % QR, &
                  loc_LHS_rows, loc_LHS_cols, &
                  name="+QR factorized Left-Hand Side")
         endif
         if (.not. allocated(settings % factorizations % QR_tau) ) then
            call f_allocate(settings % factorizations % QR_tau, &
                  loc_LHS_cols, name="Scalar factors for QR factorization")
         endif

         call parallel_factorization_QR( &
               loc_LHS, sc_desc_LHS, &
               settings % factorizations % QR, &
               settings % factorizations % QR_tau )

         call get_times(cpu_time, clock_time, unsynced=.true.)
         call output_times(deffmt, 'QR factorization', cpu_time, clock_time, defprio)

      else
         svd_on_qr = .false.

      endif ! need_qr

      if ( need_svd ) then
         call get_timestamps(cpu_time, clock_time)
         if (.not. allocated(settings % factorizations % S) ) then
            call f_allocate(settings % factorizations % S, &
                  glo_LHS_cols, name="SVD: Singular Values")
         endif
         sc_desc_VT = sc_desc_red_LHS
         if (.not. allocated(settings % factorizations % VT) ) then
            call f_allocate(settings % factorizations % VT, &
                  loc_red_LHS_rows, loc_red_LHS_cols, &
                  name="+SVD: Right Singular Vectors")
         endif
         if ( svd_on_qr ) then
            sc_desc_U = sc_desc_red_LHS
            if (.not. allocated(settings % factorizations % U) ) then
               call f_allocate(settings % factorizations % U, &
                     loc_red_LHS_rows, loc_red_LHS_cols, &
                     name="+SVD: Left Singular Vectors")
            endif
            ! copy R from QR decomposition into LHS
            sc_desc_red_LHS_in_full = sc_desc_LHS
            sc_desc_red_LHS_in_full(M_) = sc_desc_red_LHS(M_)
            call PDLASET("Lower", glo_LHS_cols-1, glo_LHS_cols-1, &
                  ZERO, ZERO, &
                  loc_LHS, 2, 1, sc_desc_red_LHS_in_full )
            call PDLACPY("Upper", glo_LHS_cols, glo_LHS_cols, &
                  settings % factorizations % QR, 1, 1, sc_desc_LHS, &
                  loc_LHS, 1, 1, sc_desc_red_LHS_in_full )
   
            call parallel_factorization_SVD( &
                  loc_LHS, sc_desc_red_LHS_in_full, &
                  settings % factorizations % U, sc_desc_U, &
                  settings % factorizations % S, &
                  settings % factorizations % VT, sc_desc_VT )
         else
            ! If SVD is preformed on full LHS, shape might have changed
            ! TODO: This is just a hotfix, should actually check if number
            ! of boundary points has changed
            if (geometry_changed .and. &
                    allocated(settings % factorizations % U)) then
               call f_deallocate(settings % factorizations % U, &
                     name="+SVD: Left Singular Vectors")
            endif

            sc_desc_U = sc_desc_LHS
            if (.not. allocated(settings % factorizations % U) ) then
               call f_allocate(settings % factorizations % U, &
                     loc_LHS_rows, loc_LHS_cols, &
                     name="+SVD: Left Singular Vectors")
            endif
   
            call parallel_factorization_SVD( &
                  loc_LHS, sc_desc_LHS, &
                  settings % factorizations % U, sc_desc_U, &
                  settings % factorizations % S, &
                  settings % factorizations % VT, sc_desc_VT )
         endif ! have qr_factorization

         ! DEBUG: coordinate scaling with interface
         write(info_str, '(4X,A,I5,A)') 'Completed SVD with ', &
                 size(settings % factorizations % S) ,' singular values.'
         call localorb_info(info_str)
         do i_sv = 1, size(settings % factorizations % S), 1
             write(info_str, '(6X,A,I5,A,E10.3)') 'Singular value ', i_sv, &
                     ' :', settings % factorizations % S(i_sv)
             call localorb_info(info_str)
         enddo
   
         call get_times(cpu_time, clock_time, unsynced=.true.)
         call output_times(deffmt, 'SVD factorization', cpu_time, clock_time, defprio)
      endif ! need_svd
   endif ! perform_factorization

   if (allocated(loc_LHS)) &
      call f_deallocate(loc_LHS, name="+Left-Hand Side")

   if (.not. iterative_solver_converged) then
      write(info_str,'(4X,A)') 'direct solution'
      call localorb_info(info_str)
      ! get timestamp
      call get_timestamps(cpu_time, clock_time)

      if ( need_qr ) then
         call parallel_solver_Q( &
               settings % factorizations % QR, sc_desc_LHS, &
               settings % factorizations % QR_tau, &
               loc_RHS, sc_desc_RHS, RSS )
      else
         RSS = ZERO
      endif
      if ( need_svd ) then
         call parallel_solver_SVD( &
               settings % factorizations % U, sc_desc_U, &
               settings % factorizations % S, &
               settings % factorizations % VT, sc_desc_VT, &
               loc_RHS, sc_desc_RHS, RSS )
      else
         call parallel_solver_R( &
               settings % factorizations % QR, sc_desc_LHS, &
               loc_RHS, sc_desc_RHS, RSS ) !TODO: actually do something with the RSS here
      endif

      call get_times(cpu_time, clock_time, unsynced=.true.)
      call output_times(deffmt, 'Direct solution', cpu_time, clock_time, defprio)
   endif ! .not. iterative_solver_converged

   ! STATISTICAL ANALYSIS 2
   call statistical_analysis_2( glo_LHS_rows_tot, glo_LHS_cols, &
         TSS, RSS, RMSD, R2, adjR2 )

   ! EXTRACTION
   ! extract reaction field coefficients from RHS matrix,
   ! which has been overwritten by the solution matrix
   call parallel_extract_rfc( loc_RHS, sc_desc_RHS, &
               continua, settings % sparsity_threshold )

   if (allocated(loc_RHS)) &
      call f_deallocate(loc_RHS, name="+Right-Hand Side")

   ! ANALYSIS

   ! get timestamp
   call get_timestamps(cpu_time, clock_time)

   ! STATISTICAL ANALYSIS
   call print_statistical_analysis_rfc( use_unit, &
                        TSS(1), RSS(1), RMSD(1), R2(1), adjR2(1) )

   ! SPARSITY ANALYSIS
   call print_sparsity_analysis_rfc( use_unit, continua(IND_RF)%basis )

   ! sum up time
   call get_times(cpu_time, clock_time)
   call output_times(deffmt, 'Time for fit analysis', &
          cpu_time, clock_time, defprio)

end subroutine calculate_reaction_field_coefficients_scalapack
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/get_scalapack_descriptors_for_SLE_setup
!  NAME
!    get_scalapack_descriptors_for_SLE_setup
!  SYNOPSIS

subroutine get_scalapack_descriptors_for_SLE_setup( blacs_context, &
         glo_LHS_rows, glo_LHS_cols, glo_RHS_rows, glo_RHS_cols, divisor, &
         loc_LHS_rows, loc_LHS_cols, sc_desc_LHS, &
         loc_RHS_rows, loc_RHS_cols, sc_desc_RHS)

!  PURPOSE
!    Get ScaLAPACK descriptors for distributed matrices in a system
!    of linear equations.
!    Of course: Parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: blacs_context
   integer, intent(in) :: glo_LHS_rows, glo_LHS_cols
   integer, intent(in) :: glo_RHS_rows, glo_RHS_cols
   integer, intent(in) :: divisor

   integer, intent(out) :: loc_LHS_rows, loc_LHS_cols, sc_desc_LHS(DLEN_)
   integer, intent(out) :: loc_RHS_rows, loc_RHS_cols, sc_desc_RHS(DLEN_)

!  INPUTS
!   o blacs_context -- BLACS context handle
!   o glo_LHS_rows -- number of rows globally in left-hand side of SLE
!   o glo_LHS_cols -- number of columns globally in left-hand side of SLE
!   o glo_RHS_rows -- number of rows globally in right-hand side of SLE
!   o glo_RHS_cols -- number of columns globally in right-hand side of SLE
!   o divisor -- block matrix such that block size is a multiple of divisor
!  OUTPUT
!   o loc_LHS_rows -- number of rows locally in left-hand side matrix
!   o loc_LHS_cols -- number of columns locally in left-hand side matrix
!   o sc_desc_LHS -- ScaLAPACK descriptor of LHS matrix
!   o loc_RHS_rows -- number of rows locally in right-hand side matrix
!   o loc_RHS_cols -- number of columns locally in right-hand side matrix
!   o sc_desc_RHS -- ScaLAPACK descriptor of RHS matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'get_scalapack_descriptors_for_SLE_setup'
   integer, parameter :: iprow0 = 0, ipcol0 = 0

   integer :: nprow, npcol, myprow, mypcol, nblk, info

   ! get BLACS grid information
   call BLACS_GRIDINFO( blacs_context, nprow, npcol, myprow, mypcol )

   if (npcol.ne.1) then
      call f_stop("Internal error: processor grid for matrix has to have"//&
                     " exactly one column", func)
   endif

   ! determine obtimal blocking factor
   nblk = ICEIL(glo_LHS_rows, divisor*nprow) * divisor

   ! get local dimensions
   loc_LHS_rows = max(1, NUMROC( glo_LHS_rows, nblk, myprow, iprow0, nprow ) )
   loc_RHS_rows = max(1, NUMROC( glo_RHS_rows, nblk, myprow, iprow0, nprow ) )
   loc_LHS_cols = glo_LHS_cols
   loc_RHS_cols = glo_RHS_cols

   ! create ScaLAPACK descriptors
   call DESCINIT( sc_desc_LHS, glo_LHS_rows, glo_LHS_cols, nblk, glo_LHS_cols, &
                  iprow0, ipcol0, blacs_context, loc_LHS_rows, info )
   call DESCINIT( sc_desc_RHS, glo_RHS_rows, glo_RHS_cols, nblk, glo_RHS_cols, &
                  iprow0, ipcol0, blacs_context, loc_RHS_rows, info )

end subroutine get_scalapack_descriptors_for_SLE_setup
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/get_scalapack_descriptors_for_SLE
!  NAME
!    get_scalapack_descriptors_for_SLE
!  SYNOPSIS

subroutine get_scalapack_descriptors_for_SLE( blacs_context, &
                  glo_LHS_rows, glo_LHS_cols, glo_RHS_rows, glo_RHS_cols, &
                  loc_LHS_rows, loc_LHS_cols, sc_desc_LHS, &
                  loc_RHS_rows, loc_RHS_cols, sc_desc_RHS )

!  PURPOSE
!    Get ScaLAPACK descriptors for distributed matrices in a system
!    of linear equations.
!    Of course: Parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: blacs_context
   integer, intent(in) :: glo_LHS_rows, glo_LHS_cols
   integer, intent(in) :: glo_RHS_rows, glo_RHS_cols

   integer, intent(out) :: loc_LHS_rows, loc_LHS_cols, sc_desc_LHS(DLEN_)
   integer, intent(out) :: loc_RHS_rows, loc_RHS_cols, sc_desc_RHS(DLEN_)
!  INPUTS
!   o blacs_context -- BLACS context handle
!   o glo_LHS_rows -- number of rows globally in left-hand side of SLE
!   o glo_LHS_cols -- number of columns globally in left-hand side of SLE
!   o glo_RHS_rows -- number of rows globally in right-hand side of SLE
!   o glo_RHS_cols -- number of columns globally in right-hand side of SLE
!  OUTPUT
!   o loc_LHS_rows -- number of rows locally in left-hand side matrix
!   o loc_LHS_cols -- number of columns locally in left-hand side matrix
!   o sc_desc_LHS -- ScaLAPACK descriptor of LHS matrix
!   o loc_RHS_rows -- number of rows locally in right-hand side matrix
!   o loc_RHS_cols -- number of columns locally in right-hand side matrix
!   o sc_desc_RHS -- ScaLAPACK descriptor of RHS matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer, parameter :: iprow0 = 0, ipcol0 = 0
   integer :: nprow, npcol, myprow, mypcol, nblk, info

   ! get BLACS grid information
   call BLACS_GRIDINFO( blacs_context, nprow, npcol, myprow, mypcol )

   ! determine obtimal blocking factor
   nblk = PILAENV(blacs_context, 'Double-precision')

   ! get local dimensions
   loc_LHS_rows = max(1, NUMROC( glo_LHS_rows, nblk, myprow, iprow0, nprow ) )
   loc_RHS_rows = max(1, NUMROC( glo_RHS_rows, nblk, myprow, iprow0, nprow ) )
   loc_LHS_cols = NUMROC( glo_LHS_cols, nblk, mypcol, ipcol0, npcol )
   loc_RHS_cols = NUMROC( glo_RHS_cols, nblk, mypcol, ipcol0, npcol )

   ! create ScaLAPACK descriptors
   call DESCINIT( sc_desc_LHS, glo_LHS_rows, glo_LHS_cols, nblk, nblk, &
                  iprow0, ipcol0, blacs_context, loc_LHS_rows, info )
   call DESCINIT( sc_desc_RHS, glo_RHS_rows, glo_RHS_cols, nblk, nblk, &
                  iprow0, ipcol0, blacs_context, loc_RHS_rows, info )

end subroutine get_scalapack_descriptors_for_SLE

!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_setup_internal_RHS_part
!  NAME
!    serial_setup_internal_RHS_part
!  SYNOPSIS

subroutine serial_setup_internal_RHS_part( continua, &
      RHS, &
      weights, &
      first_global_row )

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Serial version!
!
!    This particular subroutine enforces charge neutrality of the potential
!    in one continuum
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricContinuum), intent(in) :: continua(:)
   real(dp), intent(in) :: weights(:)
   integer, intent(in) :: first_global_row

   real(dp), intent(inout) :: RHS(:,:)

!  INPUTS
!   o continua -- all dielectric continua in the system TODO: maybe do this
!                 element-wise for all continua
!   o weights -- which portion each continuum cuts out of a sphere in lim r->infty
!   o first_global_row -- first row in the global LHS matrix that is to be
!                         filled with this block
!
!  OUTPUTS
!   o RHS -- parts of the RHS matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2019).
!  SOURCE

   integer :: n_continua, i_continuum, i_row, RHS_shape(2)
   real(dp) :: relative_charge

   i_row = first_global_row

   n_continua = size(continua)
   if (size(weights) .ne. n_continua) &
           call f_stop('Internal Error: Number of weights does not match number of continua')
   if (any(weights .lt. 0.e0_dp)) &
           call f_stop('Internal Error: Weights need to be non-negative')
   ! TODO: tolerance is quite arbitrary
   if (sum(weights) .gt. 1.001e0_dp .or. &
       sum(weights) .lt. 0.999e0_dp ) &
           call f_stop('Internal Error: Weights have to sum up to 1')
   RHS_shape = shape(RHS)
   if (RHS_shape(2) .ne. 1) &
           call f_stop('Internal Error: Unknown RHS shape')

   relative_charge = 1.e0_dp / (sum(weights*continua%eps))

   do i_continuum = 1, n_continua, 1
      select type(b => continua(i_continuum)%basis)
         class is (SolHarmBasis)
            select case(b%solharm_type)
               case(MPE_CONST%BASIS_REG)
                  ! just check if weight is 0, aferwards continue
                  if (weights(i_continuum) .ne. 0.e0_dp) &
                     call f_stop('Internal Error: Weights of bounded continua need to be 0')
               case(MPE_CONST%BASIS_IRR)

                  ! The equation to be enforced here is
                  !
                  !    sum_J c_J^(0,0) = Q [(sum_j w_j eps_j)^-1 - eps_i^-1]
                  !
                  ! with the LHS sum running over centers in continuum i,
                  ! expansion coefficients c, total charge Q, weights w and
                  ! relative permittivites eps

                  RHS(i_row,1) = charge * (relative_charge - 1.e0_dp / continua(i_continuum)%eps)
                  i_row = i_row + 1
               case default
                     call f_stop('Internal Error: Unknown solid harmonic type')
            end select ! b%solharm_type
         class default
            call f_stop('Internal Error: Enforcing charge neutrality for this basis '//&
                    'type not implemented')
      end select ! type(b => continua(i_continuum)%basis)
   enddo ! i_continuum


end subroutine serial_setup_internal_RHS_part


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_setup_RHS_part
!  NAME
!    serial_setup_RHS_part
!  SYNOPSIS

subroutine serial_setup_RHS_part( &
      interface, &
      delta_inv_eps, &
      v_hartree_callback, &
      RHS, &
      first_global_row )

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Serial version!
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(in) :: interface
   real(dp), intent(in) :: delta_inv_eps
   procedure(CBHartreePot) ::  v_hartree_callback
   integer, intent(in) :: first_global_row

   real(dp), intent(inout) :: RHS(:,:)
!  INPUTS
!   o interface -- points on the cavity surface
!   o delta_inv_eps -- difference of the inverse dielectric permittivities
!                      at the dielectric interface
!   o v_hartree_callback -- callback routine evaluating the Hartree potential
!   o first_global_row -- first row in the global RHS vector that is to be
!                         filled (starting at 1)
!  OUTPUT
!   o RHS -- part of the RHS vector
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'serial_setup_RHS_part'

   character(132) :: info_str
   integer :: info
   integer :: i_p, i_row
   integer :: n_points, n_bc
   real(dp), allocatable :: points(:,:)
   real(dp), allocatable :: v_h(:), v_h_gradient(:,:)
   logical :: need_v_hartree_gradient

   ! ASSERTIONS
   if (size(RHS,2).ne.1) then
      call f_stop("Internal error: function expects vector (second"//&
                     " dimension has size 1)", func)
   endif

   ! DIMENSIONS
   n_points = size(interface%p)
   n_bc = interface%n_bc

   call f_allocate(points, &
         1,3, lbound(interface % p, 1),ubound(interface % p, 1), &
         name="+Temporary array of points for Hartree potential evaluation")

   ! convert InterfacePoint(n)%coord(3) to real(3,n)
   call InterfacePoint_vector_extract(interface % p, &
         coords=points )

   need_v_hartree_gradient = n_bc > 2

   ! allocate Hartree potential arrays
   call f_allocate(v_h, &
         lbound(interface % p, 1),ubound(interface % p, 1), &
         name="+Evaluation of Hartree potential")
   if (need_v_hartree_gradient) then
      call f_allocate(v_h_gradient, &
            1,3, lbound(interface % p, 1),ubound(interface % p, 1), &
            name="+Evaluation of Hartree potential gradient")
   endif

   if (need_v_hartree_gradient) then
      call v_hartree_callback(n_points, points, v_h, v_h_gradient)
   else
      call v_hartree_callback(n_points, points, v_h)
   endif

   ! first two boundary conditions
   do i_p = lbound(interface % p, 1), ubound(interface % p, 1)
      i_row = first_global_row + (i_p-lbound(interface % p, 1))*n_bc
      RHS(i_row,1) = delta_inv_eps * v_h(i_p)
      RHS(i_row+1,1) = 0
   enddo ! i_p

   ! 3rd and 4th boundary condition
   if (n_bc > 2) then
      do i_p = lbound(interface % p, 1), ubound(interface % p, 1)
         i_row = first_global_row + (i_p-lbound(interface % p, 1))*n_bc
         RHS(i_row+2:i_row+3,1) = delta_inv_eps * matmul( &
               v_h_gradient(:,i_p), interface % p(i_p) % tangents(:,:) )
      enddo ! i_point
   endif

   if (allocated(points)) &
      call f_deallocate(points, &
            name="+Temporary array of points for Hartree potential evaluation")
   if (allocated(v_h)) &
      call f_deallocate(v_h, &
            name="+Evaluation of Hartree potential")
   if (allocated(v_h_gradient)) &
      call f_deallocate(v_h_gradient, &
            name="+Evaluation of Hartree potential gradient")

end subroutine serial_setup_RHS_part
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_setup_local_RHS_part
!  NAME
!    parallel_setup_local_RHS_part
!  SYNOPSIS

subroutine parallel_setup_local_RHS_part( &
      interface, &
      delta_inv_eps, &
      v_hartree_callback, &
      loc_RHS, sc_desc, &
      first_global_row )

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(in) :: interface
   real(dp), intent(in) :: delta_inv_eps
   procedure(CBHartreePot) ::  v_hartree_callback
   integer, intent(in) :: sc_desc(DLEN_)
   integer, intent(in) :: first_global_row

   real(dp), intent(inout) :: loc_RHS(:,:)
!  INPUTS
!   o interface -- points on the cavity surface
!   o delta_inv_eps -- difference of the inverse dielectric permittivities
!                      at the dielectric interface
!   o v_hartree_callback -- callback routine evaluating the Hartree potential
!   o sc_desc -- ScaLAPACK descriptor of local vector loc_RHS
!   o first_global_row -- first row in the global RHS vector that is to be
!                         filled (starting at 1)
!  OUTPUT
!   o loc_RHS -- local part of the RHS vector
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'parallel_setup_local_RHS_part'

   character(132) :: info_str
   integer :: info
   integer :: n_points, n_bc
   integer :: i_point, i_row_loc, i_row_glo
   integer :: n_row_blocks, i_row_block
   integer :: loc_row_lo, loc_row_hi, glo_row_lo, glo_row_hi
   integer :: n_rows_in_block, n_rows_in_leading_block
   integer :: n_points_in_block, i_point_in_block
   integer :: max_points_per_block
   real(dp), allocatable :: points(:,:)
   real(dp), allocatable :: v_h(:), v_h_gradient(:,:)
   logical :: need_v_hartree_gradient

   integer :: n_prow, n_pcol, my_prow, my_pcol, idummy

   ! get BLACS grid information
   call BLACS_GRIDINFO( sc_desc(CTXT_), n_prow, n_pcol, my_prow, my_pcol )

   ! DIMENSIONS
   n_points = size(interface%p)
   n_bc = interface%n_bc

   ! ASSERTIONS
   if (sc_desc(N_).ne.1) then
      call f_stop("Internal error: function expects vector (second"//&
                     " global dimension has size 1)", func)
   endif
   ! blocksize modulo #(boundary conditions) = 0
   if (mod(sc_desc(MB_), n_bc).ne.0) then
      call f_stop("Internal error: row blocking factor is not a multiple"//&
                     " of this interface's number of boundary conditions", func)
   endif
   ! matrix is not blocked in columns
   if (sc_desc(NB_).ne.sc_desc(N_)) then
      call f_stop("Internal error: matrix may not be blocked in columns"//&
                     " ", func)
   endif
   ! only one processor column in context
   if (n_pcol.ne.1) then
      call f_stop("Internal error: processor grid for matrix has to have"//&
                     " exactly one column", func)
   endif

   max_points_per_block = sc_desc(MB_) / n_bc

   call f_allocate(points, 3, max_points_per_block, &
         name="+Temporary array of points for Hartree potential evaluation")

   need_v_hartree_gradient = n_bc > 2

   ! allocate Hartree potential arrays
   call f_allocate(v_h, max_points_per_block, &
         name="+Evaluation of Hartree potential")
   if (need_v_hartree_gradient) then
      call f_allocate(v_h_gradient, 3,max_points_per_block, &
            name="+Evaluation of Hartree potential gradient")
   endif

   ! get local low indices for submatrix (index of first row)
   glo_row_lo = first_global_row
   call INFOG1L( glo_row_lo, sc_desc(MB_), n_prow, my_prow, &
                 sc_desc(RSRC_), loc_row_lo, idummy )

   ! get local high indices for submatrix (index of last row +1)
   glo_row_hi = first_global_row + n_points*n_bc
   call INFOG1L( glo_row_hi, sc_desc(MB_), n_prow, my_prow, &
                 sc_desc(RSRC_), loc_row_hi, idummy )


   ! get number of row-blocks to fill
   ! total number of blocks until upper index:
   !   ICEIL(local_row_hi-1,sc_desc(MB_))
   ! full blocks below lower index (not filled):
   !   (local_row_lo-1)/sc_desc(MB_)
   n_row_blocks = ICEIL(loc_row_hi-1,sc_desc(MB_)) - &
                        (loc_row_lo-1)/sc_desc(MB_)

   ! size of leading (partial) block
   n_rows_in_leading_block = sc_desc(MB_) - mod(loc_row_lo-1,sc_desc(MB_))
   n_rows_in_leading_block = min(n_rows_in_leading_block, loc_row_hi-loc_row_lo)


   ! MAIN

   ! local row index
   i_row_loc = loc_row_lo
   ! size of first row-block
   n_rows_in_block = n_rows_in_leading_block
   n_points_in_block = n_rows_in_block / n_bc

   ! loop over row-blocks
   do i_row_block = 1, n_row_blocks

      ! get global index of first row in block
      i_row_glo = INDXL2G( i_row_loc, sc_desc(MB_), my_prow, &
                           sc_desc(RSRC_), n_prow )

      ! translate to point index
      i_point = ICEIL(i_row_glo-glo_row_lo+1, n_bc)

      ! convert InterfacePoint(n)%coord(3) to real(3,n)
      call InterfacePoint_vector_extract( &
            interface % p(i_point:i_point-1+n_points_in_block), &
            coords=points(:,1:n_points_in_block) )

      if (need_v_hartree_gradient) then
         call v_hartree_callback(n_points_in_block, points, v_h, v_h_gradient)
      else
         call v_hartree_callback(n_points_in_block, points, v_h)
      endif

      ! loop over points in block
      do i_point_in_block = 1, n_points_in_block

         ! first two boundary conditions
         loc_RHS(i_row_loc,1) = delta_inv_eps * v_h(i_point_in_block)
         loc_RHS(i_row_loc+1,1) = 0

         ! 3rd and 4th boundary condition
         if (n_bc > 2) then
            loc_RHS(i_row_loc+2:i_row_loc+3,1) = delta_inv_eps * matmul( &
               v_h_gradient(:,i_point_in_block), &
               interface % p(i_point) % tangents(:,:) )
         endif

         ! update indices for next cycle
         i_row_loc = i_row_loc + n_bc
         i_point = i_point + 1

      enddo ! i_point_in_block

      ! size of next row-block
      ! (care is taken to treat fractional blocks correctly)
      n_rows_in_block = min(sc_desc(MB_), loc_row_hi-i_row_loc)
      n_points_in_block = n_rows_in_block / n_bc

   enddo ! i_row_block

   if (allocated(points)) &
      call f_deallocate(points, &
            name="+Temporary array of points for Hartree potential evaluation")
   if (allocated(v_h)) &
      call f_deallocate(v_h, &
            name="+Evaluation of Hartree potential")
   if (allocated(v_h_gradient)) &
      call f_deallocate(v_h_gradient, &
            name="+Evaluation of Hartree potential gradient")

end subroutine parallel_setup_local_RHS_part

!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_setup_internal_LHS_part
!  NAME
!    serial_setup_internal_LHS_part
!  SYNOPSIS

subroutine serial_setup_internal_LHS_part(continuum, LHS, first_global_row, first_global_col)

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Serial version!
!
!    This particular subroutine enforces charge neutrality of the potential
!    in one continuum
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricContinuum), intent(in) :: continuum
   integer, intent(in) :: first_global_row, first_global_col

   real(dp), intent(inout) :: LHS(:,:)
!  INPUTS
!   o continuum -- dielectric continuum
!   o first_global_row -- first row in the global LHS matrix that is to be
!                         filled with this block
!   o first_global_col -- first column in the global LHS matrix that is to be
!                         filled with this block
!  OUTPUT
!   o LHS -- part of the LHS matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2019).
!  SOURCE

   integer :: i_cond_type, i_row, i_col, i_center

   character(132) :: info_str

   i_col = first_global_col
   i_row = first_global_row

   select type(b => continuum%basis)
         class is (SolHarmBasis)
            select case(b%solharm_type)
               case(MPE_CONST%BASIS_REG)
                     ! Fast forward to next continuum
                     return
                case(MPE_CONST%BASIS_IRR)

                     ! The equation to be enforced here is
                     !
                     !    sum_J c_J^(0,0) = Q [(sum_j w_j eps_j)^-1 - eps_i^-1]
                     !
                     ! with the LHS sum running over centers in continuum i,
                     ! expansion coefficients c, total charge Q, weights w and
                     ! relative permittivites eps
                     !
                     ! Caveat: Coordinate scaling!

                     do i_center = 1, size(b%centers)
                        if (b%centers(i_center)%lmin .eq. 0) then

                           ! Unscaled coefficient Q^(l,m) of irregular basis is
                           !
                           !   rscale^(-l-1) * Q'^(l,m)
                           !
                           ! with the scaled coefficient Q'^(l,m).
                           ! Since we are only interested in l = 0, the following is sufficient

                           LHS(i_row,i_col) = 1.e0_dp / b%centers(i_center)%rscale
                        else
                           write(info_str, '(A)') " * Warning: External potential contains "//&
                                   "center with lmin>0."
                           call localorb_info(info_str)
                           write(info_str, '(A)') "            Apart from this being a bad idea, "//&
                                   "it might interfere with a later sanity check."
                           call localorb_info(info_str)
                           write(info_str, '(A)') "            If all centers have lmin>0, aims will "//&
                                   "abort with a warning about wrong matrix rows having been copied."
                           call localorb_info(info_str)
                        endif
                        i_col = i_col + SolHarmBasisCenter_size_lm(b%centers(i_center))

                     enddo ! i_center
               case default
                  call f_stop('Internal Error: Unknown solid harmonic type')
            end select ! b%solharm_type
         class default
            call f_stop('Internal Error: Enforcing charge neutrality for this basis '//&
                    'type not implemented')
   end select ! type(b => continua(i_continuum)%basis)

end subroutine serial_setup_internal_LHS_part

!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_setup_LHS_part
!  NAME
!    serial_setup_LHS_part
!  SYNOPSIS

subroutine serial_setup_LHS_part( interface, continuum, positive_side, &
               LHS, first_global_row, first_global_col )

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Serial version!
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(in) :: interface
   type(DielectricContinuum), intent(in) :: continuum
   logical, intent(in) :: positive_side
   integer, intent(in) :: first_global_row
   integer, intent(in) :: first_global_col

   real(dp), intent(inout) :: LHS(:,:)
!  INPUTS
!   o interface -- dielectric interface, holds necessary information like
!                  sampling points, boundary conditions, dielectric constants
!   o continuum -- dielectric continuum on corresponding side of the interface
!                  with the proper basis functions
!   o positive_side -- do the basis functions belong to the side of the 
!                      interface where the normal vector points to?
!   o first_global_row -- first row in the global LHS matrix that is to be
!                         filled with this block
!   o first_global_col -- first column in the global LHS matrix that is to be
!                         filled with this block
!  OUTPUT
!   o LHS -- part of the LHS matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'serial_setup_LHS_part'

   character(132) :: info_str
   integer :: info, &
              max_lmax, max_count_sh, count_sh, i_lm, &
              n_centers, i_center, n_points, i_point, n_bc, &
              i_row, i_col
   real(dp), allocatable :: cartesians(:,:), rlylm(:), drlylm(:,:)
   real(dp) :: sign, relvec(3), distsq


   ! select sign and permittivity
   if (positive_side) then
      sign = ONE
   else
      sign = -ONE
   endif

   ! DIMENSIONS
   n_points = size(interface%p)
   n_bc = interface%n_bc

   select type (b => continuum%basis)
   class is (SolHarmBasis)
      n_centers = size(b%centers)
      max_lmax = maxval(b%centers(:)%lmax)
      max_count_sh = (max_lmax+1)**2

      call initialize_cartesian_ylm(max_lmax)

      ! allocate cartesians array of corresponding size
      call f_allocate(cartesians, 1,n_max_cartesian, 0,max_lmax, name="cartesians")
      call f_allocate(rlylm, max_count_sh, name="rlylm")
      call f_allocate(drlylm, 3, max_count_sh, name="drlylm")


      ! MAIN

      ! start at offset
      i_row = first_global_row

      do i_point = 1, n_points

         ! start at offset
         i_col = first_global_col

         do i_center = 1, n_centers
            ! number of basis functions for this center
            count_sh = SolHarmBasisCenter_size_lm(b%centers(i_center))

            select case (b%solharm_type)
            case (MPE_CONST%BASIS_REG)
               call fill_centerblk_reg( n_boundary_conditions=n_bc, &
                     point=interface%p(i_point), &
                     permittivity=continuum%eps, &
                     sign=sign, &
                     center=b%centers(i_center), &
                     cartesians=cartesians, rlylm=rlylm, drlylm=drlylm, &
                     matrix=LHS(i_row:i_row+n_bc-1,i_col:i_col+count_sh-1) )
            case (MPE_CONST%BASIS_IRR)
               call fill_centerblk_irr( n_boundary_conditions=n_bc, &
                     point=interface%p(i_point), &
                     center=b%centers(i_center), &
                     permittivity=continuum%eps, &
                     sign=sign, &
                     cartesians=cartesians, rlylm=rlylm, drlylm=drlylm, &
                     matrix=LHS(i_row:i_row+n_bc-1,i_col:i_col+count_sh-1) )
            case default
               write(info_str,'(A,I2,A)') '*** ERROR: basis type no. ', &
                              b%solharm_type, ' is unknown'
               call localorb_info(info_str)
               call f_stop('Encountered unknown basis type', func)
            endselect

            ! column for next center
            i_col = i_col + count_sh
         enddo ! i_center

         ! row for next cycle
         i_row = i_row + n_bc

      enddo ! i_points

      call f_deallocate(cartesians, name="cartesians")
      call f_deallocate(rlylm, name="rlylm")
      call f_deallocate(drlylm, name="drlylm")

   class default
      call f_stop("Internal Error: Setting up the coefficient matrix is "//&
                                 "not implemented for this kind of basis")
   end select

end subroutine serial_setup_LHS_part
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_setup_local_LHS_part
!  NAME
!    parallel_setup_local_LHS_part
!  SYNOPSIS

subroutine parallel_setup_local_LHS_part( interface, continuum, &
         positive_side, part, sc_desc, first_global_row, first_global_col )

!  PURPOSE
!    Set up a part of the coefficient matrix in the systems of linear
!    equations that appear both in the reaction field factor MPE and
!    the reaction field coefficient MPE version.
!    Parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(in) :: interface
   type(DielectricContinuum), intent(in) :: continuum
   logical, intent(in) :: positive_side
   integer, intent(in) :: sc_desc(DLEN_)
   integer, intent(in) :: first_global_row
   integer, intent(in) :: first_global_col

   real(dp), intent(out) :: part(:,:)
!  INPUTS
!   o interface -- dielectric interface, holds necessary information like
!                  sampling points, boundary conditions, dielectric constants
!   o continuum -- dielectric continuum on corresponding side of the interface
!                  with the proper basis functions
!   o positive_side -- do the basis functions belong to the side of the 
!                      interface where the normal vector points to?
!   o sc_desc -- ScaLAPACK descriptor of local matrix part
!   o first_global_row -- first row in the global LHS matrix that is to be
!                         filled with this block
!   o first_global_col -- first column in the global LHS matrix that is to be
!                         filled with this block
!  OUTPUT
!   o part -- local part of the LHS submatrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'parallel_setup_local_LHS_part'

   character(132) :: info_str
   integer :: info, &
              n_prow, n_pcol, my_prow, my_pcol, idummy, &
              i_row_loc, n_rows_glo, i_row_glo, i_col_loc, &
              n_row_blocks, i_row_block, &
              loc_row_lo, loc_row_hi, glo_row_lo, glo_row_hi, &
              loc_col_lo, loc_col_hi, &
              n_rows_in_block, i_row_in_block, n_rows_in_leading_block, &
              max_lmax, max_count_sh, count_sh, n_basis_sh, n_bc, &
              n_points, i_point, n_centers, i_center
   real(dp), allocatable :: cartesians(:,:), rlylm(:), drlylm(:,:)
   real(dp) :: sign, relvec(3), distsq


   ! get BLACS grid information
   call BLACS_GRIDINFO( sc_desc(CTXT_), n_prow, n_pcol, my_prow, my_pcol )

   ! DIMENSIONS
   n_points = size(interface%p)
   n_bc = interface%n_bc

   ! select sign
   if (positive_side) then
      sign = ONE
   else
      sign = -ONE
   endif

   ! ASSERTIONS
   ! blocksize modulo #(boundary conditions) = 0
   if (mod(sc_desc(MB_), n_bc).ne.0) then
      call f_stop("Internal error: row blocking factor is not a multiple"//&
                     " of this interface's number of boundary conditions", func)
   endif
   ! matrix is not blocked in columns
   if (sc_desc(NB_).ne.sc_desc(N_)) then
      call f_stop("Internal error: matrix may not be blocked in columns"//&
                     " ", func)
   endif
   ! only one processor column in context
   if (n_pcol.ne.1) then
      call f_stop("Internal error: processor grid for matrix has to have"//&
                     " exactly one column", func)
   endif


   select type (b => continuum%basis)
   class is (SolHarmBasis)
      n_basis_sh = sum(SolHarmBasisCenter_size_lm(b%centers))
      n_centers = size(b%centers)
      max_lmax = maxval(b%centers(:)%lmax)
      max_count_sh = (max_lmax+1)**2

      call initialize_cartesian_ylm(max_lmax)

      ! allocate cartesians array of corresponding size
      call f_allocate(cartesians, 1,n_max_cartesian, 0,max_lmax, name="cartesians")
      call f_allocate(rlylm, max_count_sh, name="rlylm")
      call f_allocate(drlylm, 3, max_count_sh, name="drlylm")


      ! get local low indices for submatrix (index of first row)
      glo_row_lo = first_global_row
      call INFOG1L( glo_row_lo, sc_desc(MB_), n_prow, my_prow, &
                    sc_desc(RSRC_), loc_row_lo, idummy )

      ! get local high indices for submatrix (index of last row +1)
      glo_row_hi = first_global_row + n_points*n_bc
      call INFOG1L( glo_row_hi, sc_desc(MB_), n_prow, my_prow, &
                    sc_desc(RSRC_), loc_row_hi, idummy )

      ! index of first column
      loc_col_lo = first_global_col
      ! index of last column +1
      loc_col_hi = first_global_col + n_basis_sh

      ! get number of row-blocks to fill
      ! total number of blocks until upper index:
      !   ICEIL(loc_row_hi-1,sc_desc(MB_))
      ! full blocks below lower index (not filled):
      !   (loc_row_lo-1)/sc_desc(MB_)
      n_row_blocks = ICEIL(loc_row_hi-1,sc_desc(MB_)) - &
                           (loc_row_lo-1)/sc_desc(MB_)

      ! size of leading (partial) block
      n_rows_in_leading_block = sc_desc(MB_) - mod(loc_row_lo-1,sc_desc(MB_))
      n_rows_in_leading_block = min(n_rows_in_leading_block, loc_row_hi-loc_row_lo)


      ! MAIN

      ! local row index
      i_row_loc = loc_row_lo
      ! size of first row-block
      n_rows_in_block = n_rows_in_leading_block

      ! loop over row-blocks
      do i_row_block = 1, n_row_blocks

         ! get global index of first row in block
         i_row_glo = INDXL2G( i_row_loc, sc_desc(MB_), my_prow, &
                              sc_desc(RSRC_), n_prow )

         ! translate to point index
         i_point = ICEIL(i_row_glo-glo_row_lo+1, n_bc)

         ! loop over points in block in steps of #(boundary conditions)
         do i_row_in_block = 1, n_rows_in_block, n_bc

            ! start at first column
            i_col_loc = loc_col_lo

            do i_center = 1, n_centers
               ! number of basis functions for this center
               count_sh = SolHarmBasisCenter_size_lm(b%centers(i_center))

               select case (b%solharm_type)
               case (MPE_CONST%BASIS_REG)
                  call fill_centerblk_reg( n_boundary_conditions=n_bc, &
                        point=interface%p(i_point), &
                        permittivity=continuum%eps, &
                        sign=sign, &
                        center=b%centers(i_center), &
                        cartesians=cartesians, rlylm=rlylm, drlylm=drlylm, &
                        matrix=part(i_row_loc:i_row_loc+n_bc-1,&
                                    i_col_loc:i_col_loc+count_sh-1) )
               case (MPE_CONST%BASIS_IRR)
                  call fill_centerblk_irr( n_boundary_conditions=n_bc, &
                        point=interface%p(i_point), &
                        center=b%centers(i_center), &
                        permittivity=continuum%eps, &
                        sign=sign, &
                        cartesians=cartesians, rlylm=rlylm, drlylm=drlylm, &
                        matrix=part(i_row_loc:i_row_loc+n_bc-1,&
                                    i_col_loc:i_col_loc+count_sh-1) )
               case default
                  write(info_str,'(A,I2,A)') '*** ERROR: basis type no. ', &
                              b%solharm_type, ' is unknown'
                  call localorb_info(info_str)
                  call f_stop('Encountered unknown basis type', func)
               endselect

               ! column for next center
               i_col_loc = i_col_loc + count_sh
            enddo ! i_center

            ! update local row index for next cycle
            i_row_loc = i_row_loc + n_bc

            ! update point index for next cycle
            i_point = i_point + 1

         enddo ! i_row_in_block

         ! size of next row-block
         ! (care is taken to treat fractional blocks correctly)
         n_rows_in_block = min(sc_desc(MB_), loc_row_hi-i_row_loc)

      enddo ! i_row_block

      call f_deallocate(cartesians, name="cartesians")
      call f_deallocate(rlylm, name="rlylm")
      call f_deallocate(drlylm, name="drlylm")

   class default
      call f_stop("Internal Error: Setting up the coefficient matrix is "//&
                                 "not implemented for this kind of basis")
   end select

end subroutine parallel_setup_local_LHS_part

! <INNER AREA PROCEDURES>

!-------------------------------------------------------------------------------

subroutine fill_centerblk_irr( n_boundary_conditions, point, center, &
               permittivity, sign, cartesians, &
               rlylm, drlylm, matrix )

! USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_boundary_conditions
   type(InterfacePoint), intent(in) :: point
   type(SolHarmBasisCenter), intent(in) :: center
   real(dp), intent(in) :: permittivity
   real(dp), intent(in) :: sign
   real(dp), intent(inout) :: cartesians(:,:)
   real(dp), intent(inout) :: rlylm(:)
   real(dp), intent(inout) :: drlylm(:,:)

   real(dp), intent(inout) :: matrix(:,:)
!  INPUTS
!   o n_boundary_conditions -- number of applied boundary conditions
!   o point -- sampling point with assigned coordinate system
!   o center -- coordinates of the expansion center
!   o permittivity -- dielectric permittivity, prefactor for normal direction
!   o sign -- prefactor for all terms
!   o n_max_cartesian -- leading dimension of cartesians array
!   o cartesians -- working array for cartesian_ylm subroutines
!   o rlylm -- working array for cartesian_ylm subroutines (r^l Y_l^m)
!   o drlylm -- working array for cartesian_ylm subroutines (derivatives)
!  OUTPUT
!   o matrix -- matrix to be filled
!  SOURCE

   integer :: i_l, i_m, i_lm, i_col
   real(dp) :: dist_vec(3), inv_distance, inv_dist_sq
   real(dp) :: proj_dist_vec_on_tangents(2)
   real(dp) :: proj_dist_vec_on_normal, proj_drlylm
   real(dp) :: prefactor, prefactor_chain_rule
   real(dp) :: irr_rlylm        ! irregular solid harmonic
   real(dp) :: proj_irr_drlylm  ! projection of irregular solid harmonic's
                              !   gradient on cavity surface normal

   ! A note on prefactors:
   ! The full prefactor for solid harmonics would include
   !    sqrt( 4\pi / (2l+1) )
   ! By leaving out this part here, it is automatically assigned
   ! to the coefficients in the solution of the linear equation system.
   ! This further means that we don't have to multiply it there later.

   ! get relative distance vector of point to center
   dist_vec = (point%coord - center%coord) * center%rscale

   ! get inverse distance to center needed in prefactor to solid harmonics
   inv_dist_sq = ONE / dot_product(dist_vec,dist_vec)
   inv_distance = sqrt(inv_dist_sq)

   ! first, cartesians array has to be evaluated
   call evaluate_onecenter_cartesians( &
            dist_vec, center%lmax, cartesians)

   ! then, calculate cartesian rlylm and their gradient
   call evaluate_onecenter_cartesian_gradient_terms( &
            center%lmax, cartesians, rlylm, drlylm)
   ! apply scaling factor (chain rule)
   drlylm = drlylm * center%rscale

   ! we need the projection of the distance vector onto the
   !   surface normal to calculate the projection of the
   !   irregular solid harmonics' gradients from regular
   !   solid harmonics
   proj_dist_vec_on_normal = dot_product( point%normal, dist_vec)

   ! the same holds for the tangential directions
   proj_dist_vec_on_tangents(1) = dot_product( point%tangents(:,1), dist_vec )
   proj_dist_vec_on_tangents(2) = dot_product( point%tangents(:,2), dist_vec )

   i_lm = SolHarmBasisCenter_lbound_lm(center) - 1
   i_col = 0
   do i_l = center%lmin, center%lmax

      ! calculate prefactors
      prefactor = (inv_dist_sq**i_l)*inv_distance
      prefactor_chain_rule = center%rscale * real(2*i_l+1,dp) * inv_dist_sq

      do i_m = -i_l, i_l
         i_lm = i_lm + 1
         i_col = i_col + 1

         ! FIRST BOUNDARY CONDITION:
         !  continuity of potential

         ! calculate irregular solid harmonics from rlylm:
         !   irr_rlylm = prefactor * rlylm
         irr_rlylm = prefactor * rlylm(i_lm)

         ! assign first boundary condition (including sign)
         matrix(1,i_col) = sign * irr_rlylm

         ! done?
         if (n_boundary_conditions.le.1) cycle


         ! SECOND BOUNDARY CONDITION:
         !  continuity of electric flux in normal direction

         ! the projected derivative of the irregular solid harmonic
         ! is calculated from the projected gradient of the regular
         ! one via the chain rule of derivation

         ! projection of regular harmonic's gradient
         proj_drlylm = dot_product( point%normal, drlylm(:,i_lm) )

         ! calculate projected gradient of irregular harmonic
         proj_irr_drlylm = prefactor * &
                            ( proj_drlylm - prefactor_chain_rule * &
                              rlylm(i_lm) * proj_dist_vec_on_normal )

         ! assign second boundary condition (including sign)
         matrix(2,i_col) = sign * permittivity * proj_irr_drlylm

         ! done?
         if (n_boundary_conditions.le.2) cycle


         ! THIRD BOUNDARY CONDITION:
         !  continuity of electric field in tangential direction (I)

         ! the projected derivative of the irregular solid harmonic
         ! is calculated from the projected gradient of the regular
         ! one via the chain rule of derivation

         ! projection of regular harmonic's gradient
         proj_drlylm = dot_product( point%tangents(:,1), drlylm(:,i_lm) )

         ! calculate projected gradient of irregular harmonic
         proj_irr_drlylm = prefactor * &
                      ( proj_drlylm - prefactor_chain_rule * &
                        rlylm(i_lm) * proj_dist_vec_on_tangents(1) )

         ! assign second boundary condition (including sign)
         matrix(3,i_col) = sign * proj_irr_drlylm

         ! done?
         if (n_boundary_conditions.le.3) cycle


         ! FOURTH BOUNDARY CONDITION:
         !  continuity of electric field in tangential direction (II)

         ! the projected derivative of the irregular solid harmonic
         ! is calculated from the projected gradient of the regular
         ! one via the chain rule of derivation

         ! projection of regular harmonic's gradient
         proj_drlylm = dot_product( point%tangents(:,2), drlylm(:,i_lm) )

         ! calculate projected gradient of irregular harmonic
         proj_irr_drlylm = prefactor * &
                      ( proj_drlylm - prefactor_chain_rule * &
                        rlylm(i_lm) * proj_dist_vec_on_tangents(2) )

         ! assign second boundary condition (including sign)
         matrix(4,i_col) = sign * proj_irr_drlylm

      enddo ! i_m
   enddo ! i_l

end subroutine fill_centerblk_irr
!******
!-------------------------------------------------------------------------------

subroutine fill_centerblk_reg( n_boundary_conditions, point, permittivity, &
               sign, center, cartesians, rlylm, drlylm, matrix )

! USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_boundary_conditions
   type(InterfacePoint), intent(in) :: point
   real(dp), intent(in) :: permittivity
   real(dp), intent(in) :: sign
   type(SolHarmBasisCenter), intent(in) :: center
   real(dp), intent(inout) :: cartesians(:,:)
   real(dp), intent(inout) :: rlylm(:)
   real(dp), intent(inout) :: drlylm(:,:)

   real(dp), intent(inout) :: matrix(:,:)
!  INPUTS
!   o n_boundary_conditions -- number of applied boundary conditions
!   o point -- sampling point with assigned coordinate system
!   o permittivity -- dielectric permiitivity, prefactor for normal direction
!   o sign -- prefactor for all terms
!   o center -- expansion center with coordinates
!   o n_max_cartesian -- leading dimension of cartesians array
!   o cartesians -- working array for cartesian_ylm subroutines
!   o rlylm -- working array for cartesian_ylm subroutines (r^l Y_l^m)
!   o drlylm -- working array for cartesian_ylm subroutines (derivatives)
!  OUTPUT
!   o matrix -- matrix to be filled
!  SOURCE

   integer :: i_l, i_m, i_lm, i_col
   real(dp) :: dist_vec(3)
   real(dp) :: proj_drlylm

   ! A note on prefactors:
   ! The full prefactor for solid harmonics would include
   !    sqrt( 4\pi / (2l+1) )
   ! By leaving out this part here, it is automatically assigned
   ! to the coefficients in the solution of the linear equation system.
   ! This further means that we don't have to multiply it there later.

   ! get relative distance vector of point to center
   dist_vec = (point%coord - center%coord) * center%rscale

   ! first, cartesians array has to be evaluated
   call evaluate_onecenter_cartesians( &
            dist_vec, center%lmax, cartesians)

   ! then, calculate cartesian rlylm and their gradient
   call evaluate_onecenter_cartesian_gradient_terms( &
            center%lmax, cartesians, rlylm, drlylm)
   ! apply scaling factor (chain rule)
   drlylm = drlylm * center%rscale

   i_lm = SolHarmBasisCenter_lbound_lm(center) - 1
   i_col = 0
   do i_l = center%lmin, center%lmax
      do i_m = -i_l, i_l
         i_lm = i_lm + 1
         i_col = i_col + 1

         ! FIRST BOUNDARY CONDITION:
         !  continuity of potential

         ! assign first boundary condition (including sign)
         matrix(1,i_col) = sign * rlylm(i_lm)

         ! done?
         if (n_boundary_conditions.le.1) cycle


         ! SECOND BOUNDARY CONDITION:
         !  continuity of electric flux in normal direction

         ! projection of regular harmonic's gradient is calculated
         ! from the projection of drlylm
         proj_drlylm = dot_product( point%normal, drlylm(:,i_lm) )

         ! assign second boundary condition (including sign)
         matrix(2,i_col) = sign * permittivity * proj_drlylm

         ! done?
         if (n_boundary_conditions.le.2) cycle


         ! THIRD BOUNDARY CONDITION:
         !  continuity of electric field in tangential direction (I)

         ! projection of regular harmonic's gradient is calculated
         ! from the projection of drlylm
         proj_drlylm = dot_product( point%tangents(:,1), drlylm(:,i_lm) )

         ! assign second boundary condition (including sign)
         matrix(3,i_col) = sign * proj_drlylm

         ! done?
         if (n_boundary_conditions.le.3) cycle


         ! FOURTH BOUNDARY CONDITION:
         !  continuity of electric field in tangential direction (II)

         ! projection of regular harmonic's gradient is calculated
         ! from the projection of drlylm
         proj_drlylm = dot_product( point%tangents(:,2), drlylm(:,i_lm) )

         ! assign second boundary condition (including sign)
         matrix(4,i_col) = sign * proj_drlylm

      enddo ! i_m
   enddo ! i_l

end subroutine fill_centerblk_reg
!******

! </INNER AREA PROCEDURES>

!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_factorization_SVD
!  NAME
!    serial_factorization_SVD
!  SYNOPSIS

subroutine serial_factorization_SVD(A, U, S, VT, use_driver)

!  PURPOSE
!    Perform an incomplete singular value decomposition (SVD) of matrix A into 
!
!        A(m,n) = U(m,n) * S(n,n) * V^T(n,n)
!
!    where U and V are otrhogonal matrices and S is a diagonal matrix
!    containing the singular values.
!    Dimensions of the decomposition are determined from matrix U.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(inout) :: A(:,:)
   real(dp), intent(out) :: U(:,:)
   real(dp), intent(out) :: S(:)
   real(dp), intent(out) :: VT(:,:)
   character(len=8), intent(in), optional :: use_driver
!  INPUTS
!   o A -- matrix A
!  OUTPUT
!   o A -- has been overwritten by its bidiagonalization (see DGEBRD)
!   o U -- matrix U containing the first N left singular vectors
!   o S -- diagonal elements of matrix S sorted in descending order
!   o VT -- matrix V with dimensions (N,N) containing the right singular vectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "serial_factorization_SVD"
   character(len=8), parameter :: drivers(4) = (/ &
      "SDD     ", &
      "SVD     ", &
      "SVJ     ", &
      "custom  " /)

   character(len=8) :: driver

   character(len=80) :: info_str
   integer :: info

   integer :: M, N, LDA

   integer :: lwork
   integer :: iwork_dummy(1), ir, ic
   integer, allocatable :: iwork(:)
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: E(:), tauq(:), taup(:), work(:)

   if (present(use_driver)) then
           driver = use_driver
   else
           driver = 'custom  '
   endif

   M = size(U,1)
   N = size(U,2)
   LDA = size(A,1)

   ! ASSERTIONS
   if ( M.lt.N ) then
      call f_stop("Internal error: input matrix has wrong shape", func)
   endif
   if ( (size(S,1).ne.N) &
         .or.(LDA.lt.M).or.(size(A,2).ne.N) &
         .or.(size(VT,1).ne.N).or.(size(VT,2).ne.N) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   if (driver .eq. "SDD") then

      call f_allocate(iwork, 8*N, name="iwork")
      ! query workspace
      call DGESDD("S", M, N, A, LDA, S, U, M, VT, N, work_dummy, -1, iwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DGESSD (query) returned error ", info
         call f_stop(info_str, func)
      endif
      lwork = int(work_dummy(1))
      call f_allocate(work, lwork, name="work")

      call DGESDD("S", M, N, A, LDA, S, U, M, VT, N, work, lwork, iwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DDGESSD returned error ", info
         call f_stop(info_str, func)
      endif

   elseif (driver .eq. "SVD") then

      ! query workspace
      call DGESVD("S", "S", M, N, A, LDA, S, U, M, VT, N, work_dummy, -1, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DGESVD (query) returned error ", info
         call f_stop(info_str, func)
      endif
      lwork = int(work_dummy(1))
      call f_allocate(work, lwork, name="work")

      call DGESVD("S", "S", M, N, A, LDA, S, U, M, VT, N, work, lwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DDGESVD returned error ", info
         call f_stop(info_str, func)
      endif

   elseif (driver .eq. "SVJ") then

      lwork = max(6, M+N)
      call f_allocate(work, lwork, name="work")

      call DGESVJ("G", "U", "V", M, N, A, LDA, S, 0, U, M, work, lwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DGESVJ returned error ", info
         call f_stop(info_str, func)
      endif

      if (work(1).ne.ONE) then
         write(info_str,"(A,2X,E12.5)") "scaling of singular vectors:", work(1)
         call localorb_info(info_str)
         write(info_str,"(A,I9)") "DGESVJ detected possible underflow or overflow ", info
         call f_stop(info_str, func)
      endif

      ! VT <- V^T (where V has been stored in U)
      do ic = 1, N
         do ir = 1, N
            VT(ir,ic) = U(ic,ir)
         enddo
      enddo
      ! A contains U
      call DLACPY("F", M, N, A, LDA, U, M)

   elseif (driver .eq. "custom") then

      call f_allocate(E, N-1, name="E")
      call f_allocate(tauq, N, name="tauq")
      call f_allocate(taup, N, name="taup")
      call f_allocate(iwork, 8*N, name="iwork")

      lwork = 3 * N**2 + 4 * N ! DBDSDC
      ! query workspace
      call DGEBRD(M, N, A, LDA, S, E, tauq, taup, work_dummy, -1, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DGEBRD (query) returned error ", info
         call f_stop(info_str, func)
      endif
      lwork = max(lwork, int(work_dummy(1)))
      call DORMBR("Q", "L", "N", M, N, N, A, LDA, tauq, U, M, work_dummy, -1, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DORMBR (Q,query) returned error ", info
         call f_stop(info_str, func)
      endif
      lwork = max(lwork, int(work_dummy(1)))
      call DORMBR("P", "R", "T", N, N, M, A, LDA, taup, VT, N, work_dummy, -1, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DORMBR (P,query) returned error ", info
         call f_stop(info_str, func)
      endif
      lwork = max(lwork, int(work_dummy(1)))

      call f_allocate(work, lwork, name="work")


      call DGEBRD(M, N, A, LDA, S, E, tauq, taup, work, lwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DGEBRD returned error ", info
         call f_stop(info_str, func)
      endif

      call DBDSDC("U", "I", N, S, E, U, M, VT, N, work_dummy, iwork_dummy, &
            work, iwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DBDSDC returned error ", info
         call f_stop(info_str, func)
      endif

      call DORMBR("Q", "L", "N", M, N, N, A, LDA, tauq, U, M, work, lwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DORMBR (Q) returned error ", info
         call f_stop(info_str, func)
      endif
      call DORMBR("P", "R", "T", N, N, M, A, LDA, taup, VT, N, work, lwork, info)
      if (info.ne.0) then
         write(info_str,"(A,I9)") "DORMBR (P) returned error ", info
         call f_stop(info_str, func)
      endif

   else
         write(info_str,"(A)") "Implementation error: SVD driver routine unknown"
         call f_stop(info_str, func)
   endif ! driver

   if (allocated(iwork)) &
      call f_deallocate(iwork, name="iwork")
   if (allocated(E)) &
      call f_deallocate(E, name="E")
   if (allocated(tauq)) &
      call f_deallocate(tauq, name="tauq")
   if (allocated(taup)) &
      call f_deallocate(taup, name="taup")
   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine serial_factorization_SVD
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_solver_SVD
!  NAME
!    serial_solver_SVD
!  SYNOPSIS

subroutine serial_solver_SVD(U, S, VT, b, RSS, rcond)

!  PURPOSE
!    Solve an overdetermined system of linear equations
!
!        A(m,n) * x(n,k) = b(m,k),   m > n
!
!    using the singular value decomposition (SVD) of matrix A
!
!        A(m,n) = U(m,m) * S(m,n) * V^T(n,n)
!
!    where U and V are otrhogonal matrices and S is a diagonal matrix
!    containing the singular values.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: U(:,:)
   real(dp), intent(in) :: S(:)
   real(dp), intent(in) :: VT(:,:)
   real(dp), intent(inout) :: b(:,:)
   real(dp), intent(inout) :: RSS(:)

   real(dp), intent(in), optional :: rcond
!  INPUTS
!   o U -- matrix U with dimensions (M,N)
!   o S -- diagonal elements of matrix S sorted in descending order
!   o VT -- transposed matrix V with dimensions (N,N)
!   o b -- right-hand side of the least-square problem, dimensions (M,NRHS)
!   o RSS -- residual sum of squares from serial_solver_Q
!   o rcond -- threshold for conditioning of matrix
!  OUTPUT
!   o b -- now containing the solution vector
!   o RSS - residual sum of squares, now including the loss from discarded SVs
!           and, for non-square U, from the non-injective transformation x=U^T*b
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "serial_solver_SVD"

   character(len=132) :: info_str
   integer :: info

   integer :: M, N, LDB, NRHS

   integer :: ir, rank
   real(dp) :: Smin, rcnd
   real(dp), allocatable :: xT(:,:), Ux(:,:)

   M = size(U,1)
   N = size(U,2)
   LDB = size(b,1)
   NRHS = size(b,2)

   ! ASSERTIONS
   if ( M.lt.N ) then
      call f_stop("Internal error: input matrix has wrong shape", func)
   endif
   if ( (size(VT,1).ne.N).or.(size(VT,2).ne.N).or.(size(S,1).ne.N) &
         .or.(LDB.lt.M) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   call f_allocate(xT, NRHS, N, name="xT")

   ! x <- U^T * b
   ! but instead we save x^T <- b^T * U
   call DGEMM("T", "N", NRHS, N, M, ONE, b, LDB, U, M, ZERO, xT, NRHS)

   ! Calculate RSS. In contrast to QR-factorization, which works with elementary
   ! reflectors, here it is faster to calculate b - U*x explicitly
   ! TODO: This only holds if M > 2*N
   if (M .ne. N) then
      call f_allocate(Ux, LDB, NRHS, name='Ux')
      ! Ux <- U * x (with x = (x^T)^T)
      call DGEMM("N", "T", M, NRHS, N, ONE, U, M, xT, NRHS, ZERO, Ux, LDB)
      RSS = RSS + sum((b-Ux)**2)
      call f_deallocate(Ux, name='Ux')
   endif

   ! obtain reduced rank
   rcnd = DLAMCH('Epsilon')
   if (present(rcond)) rcnd = rcond
   Smin = rcnd * S(1)
   write(info_str, '(4X,A,E10.3)') 'Discarding singular values below ', Smin
   call localorb_info(info_str)
   ir = N
   do while ( S(ir) .lt. Smin )
      ir = ir - 1
   enddo
   rank = ir
   if (rank.lt.N) then
      write(info_str,"(6X,A,I7,A,I7,A)") &
         "The coefficient matrix of the MPE equations is rank deficient (", &
         rank, "/", N, ")."
      call localorb_info(info_str)
   else
      write(info_str,"(6X,A)") &
         "The coefficient matrix of the MPE equations has full rank."
      call localorb_info(info_str)
   endif

   ! scale x with inverted sigma
   do ir = 1, rank
      call DRSCL(NRHS, S(ir), xT(1,ir), 1)
   enddo

   ! calculate RSS
   if (NRHS .ne. 1) &
      call f_stop("Internal Error: RHS must have 1 column for statistical analysis")
   RSS = RSS + sum(xT(1,rank+1:N)**2)

   ! b <- V * x
   call DGEMM("T", "T", N, NRHS, rank, ONE, VT, N, xT, NRHS, ZERO, b, LDB)

   if (allocated(xT)) &
      call f_deallocate(xT, name="xT")

end subroutine serial_solver_SVD
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_factorization_SVD
!  NAME
!    parallel_factorization_SVD
!  SYNOPSIS

subroutine parallel_factorization_SVD(A, desc_A, U, desc_U, S, VT, desc_VT)

!  PURPOSE
!    Perform an incomplete singular value decomposition (SVD) of matrix A into 
!
!        A(m,n) = U(m,n) * S(n,n) * V^T(n,n)
!
!    where U and V are otrhogonal matrices and S is a diagonal matrix
!    containing the singular values.
!    Dimensions of the decomposition are determined from matrix U.
!
!    Parallel version!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(inout) :: A(:,:)
   integer, intent(in) :: desc_A(DLEN_)
   real(dp), intent(out) :: U(:,:)
   integer, intent(in) :: desc_U(DLEN_)
   real(dp), intent(out) :: S(:)
   real(dp), intent(out) :: VT(:,:)
   integer, intent(in) :: desc_VT(DLEN_)
!  INPUTS
!   o A -- matrix A
!   o desc_A -- ScaLAPACK descriptor for matrix A
!   o desc_U -- ScaLAPACK descriptor for matrix U
!   o desc_VT -- ScaLAPACK descriptor for matrix VT
!  OUTPUT
!   o A -- matrix A has been destroyed
!   o U -- matrix U containing the first N left singular vectors
!   o S -- diagonal elements of matrix S sorted in descending order, global!
!   o VT -- matrix V with dimensions (N,N) containing the right singular vectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "parallel_factorization_SVD"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N

   integer :: lwork
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: work(:)

   integer :: context

   M = desc_U(M_)
   N = desc_U(N_)
   context = desc_U(CTXT_)

   ! ASSERTIONS
   if ( M.lt.N ) then
      call f_stop("Internal error: input matrix has wrong shape", func)
   endif
   if ( (size(S,1).ne.N) &
         .or.(desc_A(M_).ne.M).or.(desc_A(N_).ne.N) &
         .or.(desc_VT(M_).ne.N).or.(desc_VT(N_).ne.N) ) then
      call f_stop("Internal error: global shapes of matrices do not match", &
                     func)
   endif
   if ( (desc_A(CTXT_).ne.context).or.(desc_VT(CTXT_).ne.context) ) then
      call f_stop("Internal error: all matrices have to belong to the"//&
                     " same context", func)
   endif

   ! query workspace
   call PDGESVD("V", "V", M, N, A, 1, 1, desc_A, S, U, 1, 1, desc_U, &
                  VT, 1, 1, desc_VT, work_dummy, -1, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDGESVD (query) returned error ", info
      call f_stop(info_str, func)
   endif
   lwork = int(work_dummy(1))
   call f_allocate(work, lwork, name="work")

   call PDGESVD("V", "V", M, N, A, 1, 1, desc_A, S, U, 1, 1, desc_U, &
                  VT, 1, 1, desc_VT, work, lwork, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDGESVD returned error ", info
      call f_stop(info_str, func)
   endif

   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine parallel_factorization_SVD
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_solver_SVD
!  NAME
!    parallel_solver_SVD
!  SYNOPSIS

subroutine parallel_solver_SVD(U, desc_U, S, VT, desc_VT, b, desc_b, RSS, rcond)

!  PURPOSE
!    Solve an overdetermined system of linear equations
!
!        A(m,n) * x(n,k) = b(m,k),   m > n
!
!    using the singular value decomposition (SVD) of matrix A
!
!        A(m,n) = U(m,m) * S(m,n) * V^T(n,n)
!
!    where U and V are otrhogonal matrices and S is a diagonal matrix
!    containing the singular values.
!    Parallel version!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: U(:,:)
   integer, intent(in) :: desc_U(DLEN_)
   real(dp), intent(in) :: S(:)
   real(dp), intent(in) :: VT(:,:)
   integer, intent(in) :: desc_VT(DLEN_)
   real(dp), intent(inout) :: b(:,:)
   integer, intent(in) :: desc_b(DLEN_)
   real(dp), intent(inout) :: RSS(:)

   real(dp), intent(in), optional :: rcond
!  INPUTS
!   o U -- matrix U with dimensions (M,N)
!   o desc_U -- ScaLAPACK descriptor for matrix U
!   o S -- diagonal elements of matrix S sorted in descending order, global!
!   o VT -- transposed matrix V with dimensions (N,N)
!   o desc_VT -- ScaLAPACK descriptor for matrix VT
!   o b -- right-hand side of the least-square problem, dimensions (M,NRHS)
!   o desc_b -- ScaLAPACK descriptor for matrix b
!   o RSS -- residual sum of squares from parallel_solver_Q
!   o rcond -- threshold for conditioning of matrix
!  OUTPUT
!   o b -- now containing the solution vector
!   o RSS - residual sum of squares, now including the loss from discarded SVs
!           and, for non-square U, from the non-injective transformation x=U^T*b
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "parallel_solver_SVD"

   character(len=120) :: info_str
   integer :: info

   integer :: M, N, NRHS

   integer :: ir, rank
   real(dp) :: Smin, rcnd
   integer :: nprow, npcol, myprow, mypcol, context
   integer :: desc_x(DLEN_), loc_x_rows
   real(dp), allocatable :: x(:,:), Ux(:,:)
   real(dp), allocatable :: RSS_temp(:)

   M = desc_U(M_)
   N = desc_U(N_)
   NRHS = desc_b(N_)
   context = desc_U(CTXT_)

   allocate(RSS_temp(NRHS))
   RSS_temp = 0.e0_dp

   ! ASSERTIONS
   if ( (size(S,1).ne.N).or.(desc_b(M_).lt.M) &
         .or.(desc_VT(M_).ne.N).or.(desc_VT(N_).ne.N) ) then
      call f_stop("Internal error: global shapes of matrices do not match", &
                     func)
   endif
   if ( (desc_b(CTXT_).ne.context).or.(desc_VT(CTXT_).ne.context) ) then
      call f_stop("Internal error: all matrices have to belong to the"//&
                     " same context", func)
   endif

   call BLACS_GRIDINFO( context, nprow, npcol, myprow, mypcol )
   ! create temporary matrix
   loc_x_rows = max(1, NUMROC(N, desc_b(MB_), myprow, desc_b(RSRC_), nprow))
   call DESCINIT( desc_x, N, NRHS, desc_b(MB_), desc_b(NB_), &
            desc_b(RSRC_), desc_b(CSRC_), context, loc_x_rows, info )
   call f_allocate(x, loc_x_rows, size(b,2), name="x")

   ! x <- U^T * b
   call PDGEMM("T", "N", N, NRHS, M, ONE, U, 1, 1, desc_U, b, 1, 1, desc_b, &
         ZERO, x, 1, 1, desc_x)

   ! Calculate RSS. In contrast to QR-factorization, which works with elementary
   ! reflectors, here it is faster to calculate b - U*x explicitly
   ! TODO: This only holds if M > 2*N
   if (M .ne. N) then
      call f_allocate(Ux, size(b,1), size(b,2), name='Ux')
      ! Ux <- U * x
      call PDGEMM("N", "N", M, NRHS, N, ONE, U, 1, 1, desc_U, x, 1, 1, desc_x, &
              ZERO, Ux, 1, 1, desc_b)

      ! Prevent repeated calculation of this
      Ux = b - Ux

      do ir = 1, NRHS
         call PDDOT(M, RSS_temp(ir), Ux, 1, ir, desc_b, 1, &
                                     Ux, 1, ir, desc_b, 1)
      enddo
      ! only those tasks in the corresponding column will have the correct
      ! dot product value. others have value 0. Thus divide by
      ! the number of processes in one process column and allreduce.
      RSS_temp = RSS_temp / nprow
      call sync_vector(RSS_temp, NRHS)

      RSS = RSS + RSS_temp

      RSS_temp = 0.e0_dp
      call f_deallocate(Ux, name='Ux')
   endif


   ! obtain reduced rank
   rcnd = PDLAMCH(context, 'Epsilon')
   if (present(rcond)) rcnd = rcond
   Smin = rcnd * S(1)
   write(info_str, '(4X,A,E10.3)') 'Discarding singular values below ', Smin
   call localorb_info(info_str)
   ir = N
   do while ( S(ir) .lt. Smin )
      ir = ir - 1
   enddo
   rank = ir
   if (rank.lt.N) then
      write(info_str,"(6X,A,I7,A,I7,A)") &
         "The coefficient matrix of the MPE equations is rank deficient (", &
         rank, "/", N, ")."
      call localorb_info(info_str)
   else
      write(info_str,"(6X,A)") &
         "The coefficient matrix of the MPE equations has full rank."
      call localorb_info(info_str)
   endif

   ! scale x with inverted sigma
   do ir = 1, rank
      call PDRSCL(NRHS, S(ir), x, ir, 1, desc_x, N)
   enddo

   ! calculate RSS
   do ir = 1, NRHS
      call PDDOT(N - rank, RSS_temp(ir), x, rank+1, ir, desc_x, 1, &
                                    x, rank+1, ir, desc_x, 1)
   enddo
   ! only those tasks in the corresponding column will have the correct
   ! dot product value. others have value 0. Thus divide by
   ! the number of processes in one process column and allreduce.
   RSS_temp = RSS_temp / nprow
   call sync_vector(RSS_temp, NRHS)

   RSS = RSS + RSS_temp

   ! b <- V * x
   call PDGEMM("T", "N", N, NRHS, rank, ONE, VT, 1, 1, desc_VT, &
         x, 1, 1, desc_x, ZERO, b, 1, 1, desc_b)

   if (allocated(x)) &
      call f_deallocate(x, name="x")

   deallocate(RSS_temp)

end subroutine parallel_solver_SVD
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_factorization_QR
!  NAME
!    serial_factorization_QR
!  SYNOPSIS

subroutine serial_factorization_QR(A, QR, tau)

!  PURPOSE
!    Do a QR factorization of matrix A
!
!        A(m,n) = Q(m,n) * R(n,n)
!
!    where Q is an orthogonal matrix and R is upper triagonal.
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: A(:,:)

   real(dp), intent(out) :: QR(:,:)
   real(dp), intent(out) :: tau(:)
!  INPUTS
!   o A -- matrix with dimensions (M,N)
!  OUTPUT
!   o QR -- QR form of matrix A
!   o tau -- scalar factors of elementary reflectors of the QR decimposition
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "serial_factorization_QR"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N

   integer :: lwork
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: work(:)

   M = size(A,1)
   N = size(A,2)

   ! ASSERTIONS
   if ( M.lt.N ) then
      call f_stop("Internal error: input matrix has wrong shape", func)
   endif
   if ( (size(QR,1).ne.M).or.(size(QR,2).ne.N).or.(size(tau,1).ne.N) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   ! copy A into QR
   call DLACPY("F", M, N, A, M, QR, M)

   ! workspace query
   call DGEQRF(M, N, QR, M, tau, work_dummy, -1, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "DGEQRF (query) returned error ", info
      call f_stop(info_str, func)
   endif
   lwork = int(work_dummy(1))
   call f_allocate(work, lwork, name="work")

   ! in-place QR decomposition
   call DGEQRF(M, N, QR, M, tau, work, lwork, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "DGEQRF returned error ", info
      call f_stop(info_str, func)
   endif

   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine serial_factorization_QR
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_factorization_QR
!  NAME
!    parallel_factorization_QR
!  SYNOPSIS

subroutine parallel_factorization_QR(A, desc_A, QR, tau)

!  PURPOSE
!    Do a QR factorization of A
!
!        A(m,n) = Q(m,n) * R(n,n)
!
!    where Q is an orthogonal matrix and R is upper triagonal.
!    Parallel version! 
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: A(:,:)
   integer, intent(in) :: desc_A(DLEN_)
   real(dp), intent(out) :: QR(:,:)
   real(dp), intent(out) :: tau(:)
!  INPUTS
!   o A -- matrix with dimensions (M,N)
!   o desc_A -- ScaLAPACK descriptor for matrix A
!  OUTPUT
!   o QR -- matrix A in QR form
!   o tau -- scalar factors of elementary reflectors of the QR decimposition
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "parallel_factorization_QR"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N

   integer :: lwork
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: work(:)

   integer :: context

   M = desc_A(M_)
   N = desc_A(N_)
   context = desc_A(CTXT_)

   ! ASSERTIONS
   if ( M.lt.N ) then
      call f_stop("Internal error: input matrix has wrong shape", func)
   endif
   if ( (size(QR,1).ne.size(A,1)).or.(size(QR,2).ne.size(A,2)) &
         .or.(size(tau,1).ne.size(A,2)) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   ! copy A into QR
   call PDLACPY("F", M, N, A, 1, 1, desc_A, QR, 1, 1, desc_A)

   ! workspace query
   call PDGEQRF(M, N, QR, 1, 1, desc_A, tau, work_dummy, -1, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDGEQRF (query) returned error ", info
      call f_stop(info_str, func)
   endif
   lwork = int(work_dummy(1))
   call f_allocate(work, lwork, name="work")

   ! in-place QR decomposition
   call PDGEQRF(M, N, QR, 1, 1, desc_A, tau, work, lwork, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDGEQRF returned error ", info
      call f_stop(info_str, func)
   endif

   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine parallel_factorization_QR
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_solver_Q
!  NAME
!    serial_solver_Q
!  SYNOPSIS

subroutine serial_solver_Q(QR, tau, b, RSS)

!  PURPOSE
!    Shift orthogonal matrix Q from QR decomposition to right hand side, i.e.
!
!    Q * R * x = b  -->  R * x = Q^T * b
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: QR(:,:)
   real(dp), intent(in) :: tau(:)
   real(dp), intent(inout) :: b(:,:)

   real(dp), intent(out) :: RSS(:)
!  INPUTS
!   o QR -- left-hand side matrix A with dimensions (M,N) in QR form
!   o tau -- scalar factors of elementary reflectors of the QR decomposition
!   o b -- right-hand side of the least-square problem, dimensions (M,NRHS)
!  OUTPUT
!   o b -- now contains Q^T * b
!   o RSS -- residual sum of squares
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = "serial_solver_Q"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N, NRHS

   integer :: lwork
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: work(:)

   M = size(QR,1)
   N = size(QR,2)
   NRHS = size(b,2)

   ! ASSERTIONS
   if ( (size(b,1).ne.M).or.(size(tau,1).ne.N) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   ! query workspace
   call DORMQR('L', 'T', M, NRHS, N, QR, M, tau, b, M, work_dummy, -1, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "DORMQR (query) returned error ", info
      call f_stop(info_str, func)
   endif
   lwork = int(work_dummy(1))
   call f_allocate(work, lwork, name="work")

   ! b <- Q^T * b
   call DORMQR('L', 'T', M, NRHS, N, QR, M, tau, b, M, work, lwork, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "DORMQR returned error ", info
      call f_stop(info_str, func)
   endif

   ! calculate RSS
   RSS = sum(b(N+1:M,:)**2, dim=1)

   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine serial_solver_Q
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_solver_R
!  NAME
!    serial_solver_R
!  SYNOPSIS

subroutine serial_solver_R(QR, QTb, RSS)

!  PURPOSE
!    Solve an overdetermined system of linear equations.
!
!        R(n,n) * x(n,k) = Q(m,n)^T * b(m,k),   m > n
!
!    using an existing QR decomposition.
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: QR(:,:)
   real(dp), intent(inout) :: QTb(:,:)
   real(dp), intent(inout) :: RSS(:)
!  INPUTS
!   o QR -- left-hand side matrix A with dimensions (M,N) in QR form
!   o QTb -- right-hand side of the above equation, dimensions (N,NRHS)
!   o RSS -- residual sum of squares from serial_solver_Q
!  OUTPUT
!   o QTb -- now contains the solution to the least-square problem
!   o RSS -- TODO: This is a dummy! Will currently leave the subroutine unchanged
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = "serial_solver_QR"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N, NRHS

   M = size(QR,1)
   N = size(QR,2)
   NRHS = size(QTb,2)

   ! ASSERTIONS
   if ( (size(QTb,1).lt.N) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif

   ! backward substitution: solve R * x = Q^T * b
   call DTRTRS('U', 'N', 'N', N, NRHS, QR, M, QTb, M, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "DTRTRS returned error ", info
      call f_stop(info_str, func)
   endif
   ! now, QTb contains the solution vector (first N entries)

end subroutine serial_solver_R
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_solver_Q
!  NAME
!    parallel_solver_Q
!  SYNOPSIS

subroutine parallel_solver_Q(QR, desc_QR, tau, b, desc_b, RSS)

!  PURPOSE
!    Shift orthogonal matrix Q from QR decomposition to right hand side, i.e.
!
!    Q * R * x = b  -->  R * x = Q^T * b
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: QR(:,:)
   integer, intent(in) :: desc_QR(DLEN_)
   real(dp), intent(in) :: tau(:)
   real(dp), intent(inout) :: b(:,:)
   integer, intent(in) :: desc_b(DLEN_)

   real(dp), intent(out) :: RSS(:)
!  INPUTS
!   o QR -- left-hand side matrix A with dimensions (M,N) in QR form
!   o desc_QR -- ScaLAPACK descriptor for matrix A
!   o tau -- scalar factors of elementary reflectors of the QR decomposition
!   o b -- right-hand side of the least-square problem, dimensions (M,NRHS)
!   o desc_b -- ScaLAPACK descriptor for matrix b
!  OUTPUT
!   o b -- now contains Q^T * b
!   o RSS -- residual sum of squares
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = "parallel_solver_Q"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N, NRHS

   integer :: lwork
   real(dp) :: work_dummy(1)
   real(dp), allocatable :: work(:)

   integer :: i_rhs
   integer :: context, nprow, npcol, myprow, mypcol

   M = desc_QR(M_)
   N = desc_QR(N_)
   NRHS = desc_b(N_)
   context = desc_QR(CTXT_)

   ! get BLACS grid information
   call BLACS_GRIDINFO( context, nprow, npcol, myprow, mypcol )

   ! ASSERTIONS
   if ( (desc_b(M_).ne.M) ) then
      call f_stop("Internal error: global shapes of matrices do not match", &
                     func)
   endif
   if ( (size(tau,1).ne.size(QR,2)).or.(size(RSS,1).ne.NRHS) ) then
      call f_stop("Internal error: shapes of matrices do not match", func)
   endif
   if ( (desc_b(CTXT_).ne.context) ) then
      call f_stop("Internal error: all matrices have to belong to the"//&
                     " same context", func)
   endif

   ! query workspace
   call PDORMQR('L', 'T', M, NRHS, N, QR, 1, 1, desc_QR, tau, &
                  b, 1, 1, desc_b, work_dummy, -1, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDORMQR (query) returned error ", info
      call f_stop(info_str, func)
   endif
   lwork = int(work_dummy(1))
   call f_allocate(work, lwork, name="work")

   ! b <- Q^T * b
   call PDORMQR('L', 'T', M, NRHS, N, QR, 1, 1, desc_QR, tau, &
                  b, 1, 1, desc_b, work, lwork, info)
   if (info.ne.0) then
      write(info_str,"(A,I9)") "PDORMQR returned error ", info
      call f_stop(info_str, func)
   endif

   ! calculate RSS
   do i_rhs = 1, NRHS
      call PDDOT(M - N, RSS(i_rhs), b, N+1, i_rhs, desc_b, 1, &
                                    b, N+1, i_rhs, desc_b, 1)
   enddo
   ! only those tasks in the corresponding column will have the correct
   ! dot product value. others have value 0. Thus divide by
   ! the number of processes in one process column and allreduce.
   RSS = RSS / nprow
   call sync_vector(RSS, NRHS)

   if (allocated(work)) &
      call f_deallocate(work, name="work")

end subroutine parallel_solver_Q
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_solver_R
!  NAME
!    parallel_solver_R
!  SYNOPSIS

subroutine parallel_solver_R(QR, desc_QR, QTb, desc_QTb, RSS)

!  PURPOSE
!    Solve an overdetermined system of linear equations.
!
!        R(n,n) * x(n,k) = Q(m,n)^T * b(m,k),   m > n
!
!    using an existing QR decomposition.
!    Parallel version!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: QR(:,:)
   integer, intent(in) :: desc_QR(DLEN_)
   real(dp), intent(inout) :: QTb(:,:)
   integer, intent(in) :: desc_QTb(DLEN_)
   real(dp), intent(inout) :: RSS(:)
!  INPUTS
!   o QR -- left-hand side matrix A with dimensions (M,N) in QR form
!   o desc_QR -- ScaLAPACK descriptor for matrix A
!   o QTb -- right-hand side of the above equation, dimensions (N,NRHS)
!   o desc_QTb -- ScaLAPACK descriptor for matrix QTb
!   o RSS -- residual sum of squares from parallel_solver_Q
!  OUTPUT
!   o QTb -- now containing the solution vector
!   o RSS -- TODO: This is a dummy! Will currently leave the subroutine unchanged
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = "parallel_solver_R"

   character(len=80) :: info_str
   integer :: info

   integer :: M, N, NRHS

   integer :: context

   M = desc_QR(M_)
   N = desc_QR(N_)
   NRHS = desc_QTb(N_)
   context = desc_QR(CTXT_)

   ! ASSERTIONS
   if ( (desc_QTb(M_).lt.N) ) then
      call f_stop("Internal error: global shapes of matrices do not match", &
                     func)
   endif
   if ( (desc_QTb(CTXT_).ne.context) ) then
      call f_stop("Internal error: all matrices have to belong to the"//&
                     " same context", func)
   endif

   ! backward substitution: solve R * x = Q^T * b
   call PDTRSM('L', 'U', 'N', 'N', N, NRHS, 1.e0_dp, QR, 1, 1, desc_QR, &
                  QTb, 1, 1, desc_QTb)
   ! now, QTb contains the solution vector

end subroutine parallel_solver_R
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_statistical_analysis_1
!  NAME
!    serial_statistical_analysis_1
!  SYNOPSIS

subroutine serial_statistical_analysis_1( RHS, TSS )

!  PURPOSE
!    Statistical analysis, part I, of the reaction field coefficients and
!    reaction field factor calculation routines.
!    Serial version!
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: RHS(:,:)

   real(dp), intent(out) :: TSS(:)
!  INPUTS
!   o RHS -- RHS matrix
!   o n_rows -- rows of RHS matrix
!  OUTPUT
!   o TSS -- variation of columns
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'serial_statistical_analysis_1'

   integer :: n_rows, n_cols, i_col
   real(dp), allocatable :: mean(:)

   n_rows = size(RHS,1)
   n_cols = size(RHS,2)

   ! ASSERTIONS
   if (n_cols.ne.size(TSS,1)) then
      call f_stop("Internal error: shape mismatch", func)
   endif

   call f_allocate(mean, n_cols, name="mean")

   ! calculate mean value
   mean = sum(RHS, dim=1) / n_rows

   ! calculate total sum of squares
   do i_col = 1, n_cols
      TSS(i_col) = sum( (RHS(:,i_col) - mean(i_col))**2 )
   enddo

   call f_deallocate(mean, name="mean")

end subroutine serial_statistical_analysis_1
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_statistical_analysis_1
!  NAME
!    parallel_statistical_analysis_1
!  SYNOPSIS

subroutine parallel_statistical_analysis_1(loc_RHS, sc_desc, TSS)

!  PURPOSE
!    Statistical analysis, part I, of the reaction field coefficients
!    calculation.
!    Parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: sc_desc(DLEN_)
   real(dp), intent(in) :: loc_RHS(:,:)

   real(dp), intent(out) :: TSS(:)
!  INPUTS
!   o loc_RHS -- local part of RHS matrix
!   o sc_desc -- ScaLAPACK descriptor of RHS matrix
!  OUTPUT
!   o TSS -- variation of columns
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'parallel_statistical_analysis_1'

   integer :: n_rows, i_row, n_cols, i_col, nprow, npcol, myprow, mypcol
   real(dp), allocatable :: mean(:), temp(:,:)

   n_rows = sc_desc(M_)
   n_cols = sc_desc(N_)

   ! ASSERTIONS
   if (n_cols.ne.size(TSS,1)) then
      call f_stop("Internal error: shape mismatch", func)
   endif

   call f_allocate(mean, n_cols, name="mean")
   call f_allocate(temp, size(loc_RHS,1), size(loc_RHS,2), name="temp")

   ! get BLACS grid information
   call BLACS_GRIDINFO(sc_desc(CTXT_), nprow, npcol, myprow, mypcol)

   ! The following part is rather ugly, but unfortunately, 
   ! I couldn't find any parallel sum, so everything is done with the
   ! dot product.

   ! calculate mean value
   temp = sqrt(abs(loc_RHS))
   do i_col = 1, n_cols
      call PDDOT(n_rows, mean(i_col), temp, 1, i_col, sc_desc, 1, &
                        sign(temp,loc_RHS), 1, i_col, sc_desc, 1)
   enddo
   mean = mean / (nprow*n_rows)
   call sync_vector(mean, n_cols)

   ! calculate total sum of squares
   do i_col = 1, n_cols
      temp = loc_RHS - mean(i_col)
      call PDDOT(n_rows, TSS(i_col), temp, 1, i_col, sc_desc, 1, &
                                     temp, 1, i_col, sc_desc, 1)
   enddo
   TSS = TSS / nprow
   call sync_vector(TSS, n_cols)

   call f_deallocate(mean, name="mean")
   call f_deallocate(temp, name="temp")

end subroutine parallel_statistical_analysis_1
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/statistical_analysis_2
!  NAME
!    statistical_analysis_2
!  SYNOPSIS

subroutine statistical_analysis_2( &
      n_rows, n_variables, TSS, RSS, &
      RMSD, R2, adjR2 )

!  PURPOSE
!    Statistical analysis, part II, of the reaction field coefficients and
!    reaction field factor calculation routines.
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_variables
   real(dp), intent(in) :: TSS(:)
   real(dp), intent(in) :: RSS(:)

   real(dp), intent(out) :: RMSD(:)
   real(dp), intent(out) :: R2(:)
   real(dp), intent(out) :: adjR2(:)
!  INPUTS
!   o n_rows -- number of rows in RHS matrix
!   o n_variables -- number of variables in model
!   o TSS -- variation of columns
!   o RSS -- the already calculated sum of squared residuals
!  OUTPUT
!   o RMSD -- root mean squared deviation
!   o R2 -- coefficient of determination
!   o adjR2 -- adjusted coefficient of determination
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'statistical_analysis_2'

   ! calculate root mean squared deviation
   RMSD = sqrt( RSS/n_rows )

   ! calculate coefficient of determination R^2
   R2 = 1.e0_dp - RSS/TSS

   ! calculate adjusted coefficient of determination
   adjR2 = 1.e0_dp - ( (n_rows-1)*RSS ) / &
                         ( (n_rows-n_variables)*TSS )

end subroutine statistical_analysis_2
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/print_statistical_analysis_rfc
!  NAME
!    print_statistical_analysis_rfc
!  SYNOPSIS

subroutine print_statistical_analysis_rfc( unit,  &
                           TSS, RSS, RMSD, R2, adjR2 )

!  PURPOSE
!    Print a short (or long) resume of the statistical analysis
!    of the reaction field coefficient calculation.
!    Serial AND parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: unit
   real(dp), intent(in) :: TSS, RSS, RMSD, R2, adjR2

!  INPUTS
!   o unit -- unit for output
!   o RSS -- sum of squared residuals
!   o RMSD -- root mean squared deviation
!   o R2 -- coefficient of determination
!   o adjR2 -- adjusted coefficient of determination
!  OUTPUT
!    prints to specified output
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(132) :: info_str

   ! FIT QUALITY
   write(info_str,'(4X,A)') 'Statistical analysis:'
   call localorb_info(info_str)
   write(info_str,'(8X,A30,1X,E12.5)') 'RMSD', RMSD
   call localorb_info(info_str)
   write(info_str,'(8X,A30,1X,F12.4)') 'R2', R2
   call localorb_info(info_str)
   write(info_str,'(8X,A30,1X,F12.4)') 'adjR2', adjR2
   call localorb_info(info_str)

end subroutine print_statistical_analysis_rfc
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/print_sparsity_analysis_rfc
!  NAME
!    print_sparsity_analysis_rfc
!  SYNOPSIS

subroutine print_sparsity_analysis_rfc( unit, b )

!  PURPOSE
!    Print a short resume of the sparsity analysis
!    of the reaction field coefficient calculation.
!    Serial AND parallel version!
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: unit
   class(Basis), intent(in) :: b

!  INPUTS
!   o unit -- unit for output
!   o b -- includes centers with expansion coefficients
!  OUTPUT
!    prints to specified output
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(132) :: info_str
   integer :: n_centers, i_center, n_elements, n_nonzero


   ! SPARSITY
   write(info_str,'(4X,A)') 'Sparsity analysis'
   call localorb_info(info_str)

   ! overview
   n_elements = b%get_size()

   select type (b)
   class is (SolHarmBasis)
      n_centers = size(b%centers)

      ! count number of non-zero elements
      n_nonzero = 0
      do i_center = 1, n_centers
         n_nonzero = n_nonzero + size(b%centers(i_center)%coeff)
      enddo ! i_center

      write(info_str,'(8X,A30,1X,I12)') &
            'total number of elements', n_elements
      call localorb_info(info_str)
      write(info_str,'(8X,A30,1X,I12)') &
            'number of non-zero elements', n_nonzero
      call localorb_info(info_str)
      write(info_str,'(8X,A30,1X,E12.5)') &
            'over all sparsity factor', n_nonzero/real(n_elements,dp)
      call localorb_info(info_str)

      ! broken down by expansion center
      write(info_str,'(6X,A)') 'by expansion centers'
      call localorb_info(info_str)
      write(info_str,'(8X,A6,1X,A8,1X,A8,1X,A10)') 'center', &
                        'elements', 'non-zero', 'factor'
      call localorb_info(info_str)

      do i_center = 1, n_centers
         n_elements = SolHarmBasisCenter_size_lm(b%centers(i_center))
         n_nonzero = size(b%centers(i_center)%coeff)
         write(info_str,'(8X,I6,1X,I8,1X,I8,1X,E10.3)') i_center, n_elements, &
                                 n_nonzero, n_nonzero/real(n_elements,dp)
         call localorb_info(info_str)
      enddo ! i_center

   class default
      write(info_str,'(6X,A)') "is either not possible or not implemented "//&
                              "for this kind of basis"
      call localorb_info(info_str)
   end select

end subroutine print_sparsity_analysis_rfc
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/extract_single_center_solution_reg
!  NAME
!    extract_single_center_solution_reg
!  SYNOPSIS

subroutine extract_single_center_solution_reg( dense, center, threshold )

!  PURPOSE
!    Extract the basis coefficients from a part of the solution vector and
!    store the in a compressed sparse format.
!
!    ANNOTATION:
!    The reaction field factors' prefactor comprises of
!    different contributions:
!     - the multipole moments in aims are scaled by a factor
!       of 1/sqrt(4pi), thus we need to multiply with
!       sqrt(4pi)
!     - the regular spherical harmonics carry a prefactor
!       of sqrt(4pi/(2l'+1)) that is already multiplied to
!       the reaction field coefficients to save computations
!     - the fact that the interaction energy is only 1/2 of
!       the product of "M * f_llmm * M" is accounted for by
!       correcting the total energy in the end
!   
!    This prefactor has been verified by
!     1. comparing the results of atomic ions with the analytical
!        Born equations (constant term)
!     2. checking the consistency of results of single centered
!        vs. multicentered expansions of the reaction field
!        (l-dependence)
!   
!    However, we don't need to do any scaling here because the
!    correct prefactors are implicitly assigned by definition
!    of the coefficients in the left-hand side of the equation system.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: dense(:)
   type(SolHarmBasisCenter), intent(inout) :: center
   real(dp), intent(in) :: threshold
!  INPUTS
!   o dense -- part of the solution vector belonging to center
!   o center -- the basis center
!   o threshold -- threshold for a matrix entry to be considered zero
!  OUTPUT
!   o center -- the basis center with expansion coefficients
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: i_l, i_m, i_dense
   real(dp), allocatable :: scaled_solution(:)
   real(dp) :: rscale, scaling

   rscale = center%rscale
   scaling = ONE
   do i_l = 0, center%lmin-1
      scaling = scaling * rscale
   enddo ! i_l

   call f_allocate(scaled_solution, size(dense), name="scaled_solution")
   i_dense = 0
   do i_l = center%lmin, center%lmax
      do i_m = -i_l, i_l
         i_dense = i_dense + 1
         scaled_solution(i_dense) = dense(i_dense)*scaling
      enddo ! i_m
      scaling = scaling * rscale
   enddo ! i_l

   call SpVecTuple_vector_from_dense_vector( scaled_solution, &
                  SolHarmBasisCenter_lbound_lm(center)-1, &
                  threshold, center%coeff )

   call f_deallocate(scaled_solution, name="scaled_solution")

end subroutine extract_single_center_solution_reg
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/extract_single_center_solution_irr
!  NAME
!    extract_single_center_solution_irr
!  SYNOPSIS

subroutine extract_single_center_solution_irr( dense, center, threshold )

!  PURPOSE
!    Extract the basis coefficients from a part of the solution vector and
!    store the in a compressed sparse format.
!
!    ANNOTATION:
!    The reaction field factors' prefactor comprises of
!    different contributions:
!     - the multipole moments in aims are scaled by a factor
!       of 1/sqrt(4pi), thus we need to multiply with
!       sqrt(4pi)
!     - the regular spherical harmonics carry a prefactor
!       of sqrt(4pi/(2l+1)) that is already multiplied to
!       the reaction field coefficients to save computations
!     - the fact that the interaction energy is only 1/2 of
!       the product of "M * f_llmm * M" is accounted for by
!       correcting the total energy in the end
!   
!    This prefactor has been verified by
!     1. comparing the results of atomic ions with the analytical
!        Born equations (constant term)
!     2. checking the consistency of results of single centered
!        vs. multicentered expansions of the reaction field
!        (l-dependence)
!   
!    However, we don't need to do any scaling here because the
!    correct prefactors are implicitly assigned by definition
!    of the coefficients in the left-hand side of the equation system.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: dense(:)
   type(SolHarmBasisCenter), intent(inout) :: center
   real(dp), intent(in) :: threshold
!  INPUTS
!   o dense -- part of the solution vector belonging to center
!   o center -- the basis center
!   o threshold -- threshold for a matrix entry to be considered zero
!  OUTPUT
!   o center -- the basis center with expansion coefficients
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: i_l, i_m, i_dense
   real(dp), allocatable :: scaled_solution(:)
   real(dp) :: inv_rscale, scaling

   inv_rscale = ONE/center%rscale
   scaling = inv_rscale
   do i_l = 0, center%lmin-1
      scaling = scaling * inv_rscale
   enddo ! i_l

   call f_allocate(scaled_solution, size(dense), name="scaled_solution")
   i_dense = 0
   do i_l = center%lmin, center%lmax
      do i_m = -i_l, i_l
         i_dense = i_dense + 1
         scaled_solution(i_dense) = dense(i_dense)*scaling
      enddo ! i_m
      scaling = scaling * inv_rscale
   enddo ! i_l

   call SpVecTuple_vector_from_dense_vector( scaled_solution, &
                  SolHarmBasisCenter_lbound_lm(center)-1, &
                  threshold, center%coeff )

   call f_deallocate(scaled_solution, name="scaled_solution")

end subroutine extract_single_center_solution_irr
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/serial_extract_rfc
!  NAME
!    serial_extract_rfc
!  SYNOPSIS

subroutine serial_extract_rfc( RHS, continua, threshold )

!  PURPOSE
!    Extract the reaction field coefficients from the solution vector and
!    store the in a compressed sparse format.
!    Serial version!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: RHS(:,:)
   type(DielectricContinuum), intent(inout) :: continua(:)
   real(dp), intent(in) :: threshold
!  INPUTS
!   o RHS -- full solution vector
!   o continua -- dielectric continua with basis
!   o threshold -- threshold for a matrix entry to be considered zero
!  OUTPUT
!   o continua -- dielectric continua with basis and coefficients
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'serial_extract_rfc'

   integer :: i_row, i_dc, i_center, n_elements

   character(132) :: info_str

   real(dp) :: total_coeff

   ! ASSERTIONS
   if (size(RHS,2).ne.1) then
      call f_stop("Internal error: function expects vector (second"//&
                     " dimension has size 1)", func)
   endif

   i_row = 1
   do i_dc = lbound(continua,1), ubound(continua,1)
      select type (b => continua(i_dc)%basis)
      class is (SolHarmBasis)
         total_coeff = 0.e0_dp
         do i_center = lbound(b%centers,1), ubound(b%centers,1)
            n_elements = SolHarmBasisCenter_size_lm( &
                              b%centers(i_center) )
            select case(b%solharm_type)
               case(MPE_CONST%BASIS_REG)
                  call extract_single_center_solution_reg( &
                              RHS(i_row:i_row-1+n_elements,1), &
                              b%centers(i_center), threshold )
               case(MPE_CONST%BASIS_IRR)
                  call extract_single_center_solution_irr( &
                              RHS(i_row:i_row-1+n_elements,1), &
                              b%centers(i_center), threshold )
                  total_coeff = total_coeff + b%centers(i_center)%coeff(1)%val

            end select
            i_row = i_row + n_elements
         enddo ! i_center

      class default
         call f_stop("Internal error: extraction of solution not "//&
                              "implemented for given basis type.", func)
      end select
   enddo ! i_dc

end subroutine serial_extract_rfc
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/parallel_extract_rfc
!  NAME
!    parallel_extract_rfc
!  SYNOPSIS

subroutine parallel_extract_rfc( loc_RHS, sc_desc, continua, threshold )

!  PURPOSE
!    Extract the reaction field coefficients from the solution vector and
!    store the in a compressed sparse format.
!    Parallel version!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: sc_desc(DLEN_)
   real(dp), intent(in) :: loc_RHS(:,:)
   type(DielectricContinuum), intent(inout) :: continua(:)
   real(dp), intent(in) :: threshold
!  INPUTS
!   o loc_RHS -- local part of RHS vector
!   o sc_desc -- ScaLAPACK descriptor of RHS vector
!   o continua -- dielectric continua with basis
!   o threshold -- threshold for a matrix entry to be considered zero
!  OUTPUT
!   o continua -- dielectric continua with basis and coefficients
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'parallel_extract_rfc'

   integer :: myprow, mypcol, nprow, npcol, dummy, &
              i_center, i_dc, n_elements, &
              start_row_loc, next_row_loc, i_row_loc, i_row_glo, i_row_glo_rel
   real(dp), allocatable :: dense(:)

   character(132) :: info_str
   integer :: i_coeff
   real(dp) :: total_coeff

   ! ASSERTIONS
   if (sc_desc(N_).ne.1) then
      call f_stop("Internal error: function expects vector (second"//&
                     " dimension has size 1)", func)
   endif

   ! get info about BLACS grid
   call BLACS_GRIDINFO(sc_desc(CTXT_), nprow, npcol, myprow, mypcol)


   ! get global starting row for reaction field coefficients
   i_row_glo = 1
   ! translate to local starting row
   call INFOG1L(i_row_glo, sc_desc(MB_), &
                  nprow, myprow, sc_desc(RSRC_), next_row_loc, dummy)

   do i_dc = lbound(continua,1), ubound(continua,1)
      select type (b => continua(i_dc)%basis)
      class is (SolHarmBasis)
         call f_allocate(dense, &
               maxval(SolHarmBasisCenter_size_lm(b%centers)), &
               name="dense")

         ! loop over reaction field centers
         total_coeff = 0.e0_dp
         do i_center = lbound(b%centers,1),ubound(b%centers,1)

            ! number of basis functions for center
            n_elements = SolHarmBasisCenter_size_lm(b%centers(i_center))

            ! get local start and end
            start_row_loc = next_row_loc
            call INFOG1L(i_row_glo + n_elements, sc_desc(MB_), &
                           nprow, myprow, sc_desc(RSRC_), next_row_loc, dummy)

            dense = ZERO
            ! Is my process on first column? If so, copy data
            if (mypcol.eq.sc_desc(CSRC_)) then
               do i_row_loc = start_row_loc, next_row_loc-1
                  i_row_glo_rel = INDXL2G(i_row_loc, sc_desc(MB_), myprow, &
                                          sc_desc(RSRC_), nprow) - i_row_glo + 1
                  dense(i_row_glo_rel) = loc_RHS(i_row_loc,1)
               enddo
            endif ! process is on first column

            ! synchronize dense row
            call sync_vector(dense, n_elements)

            select case(b%solharm_type)
               case(MPE_CONST%BASIS_REG)
                  call extract_single_center_solution_reg( dense(1:n_elements), &
                              b%centers(i_center), threshold )
!                  write(info_str,'(4X,A,I3,A,I1)') 'Coefficient vector at regular center ', i_center, &
!                                              ' of basis ', i_dc
!                  call localorb_info(info_str)
!                  do i_coeff = 1, size(b%centers(i_center)%coeff)
!                     write(info_str, '(6X,E15.8)') b%centers(i_center)%coeff(i_coeff)%val
!                     call localorb_info(info_str)
!                  enddo

               case(MPE_CONST%BASIS_IRR)
                  call extract_single_center_solution_irr( dense(1:n_elements), &
                              b%centers(i_center), threshold )
                  write(info_str,'(4X,A,I3,A,I1,A,E15.8)') 'Coefficient 0 at irregular center ', i_center, &
                                              ' of basis ', i_dc, ': ', b%centers(i_center)%coeff(1)%val
                  call localorb_info(info_str)
!                  do i_coeff = 1, size(b%centers(i_center)%coeff)
!                     write(info_str, '(6X,E15.8)') b%centers(i_center)%coeff(i_coeff)%val
!                     call localorb_info(info_str)
!                  enddo
                  total_coeff = total_coeff + b%centers(i_center)%coeff(1)%val
            endselect

            ! shift i_row_glo to next block
            i_row_glo = i_row_glo + n_elements

         enddo ! i_center
         write(info_str, '(2X,A,E15.8)') 'Total: ', total_coeff
         call localorb_info(info_str)

         call f_deallocate(dense, name="dense")

      class default
         call f_stop("Internal error: extraction of solution not "//&
                              "implemented for given basis type.", func)
      end select
   enddo ! i_dc

end subroutine parallel_extract_rfc
!******
!!-------------------------------------------------------------------------------
!!****s* mpe_reaction_field/sort_centers
!!  NAME
!!    sort_centers
!!  SYNOPSIS
!
!subroutine sort_centers(point, expansion_centers, &
!                        distsqs, relvecs, sorted_indices)
!
!!  PURPOSE
!!    Calculate the inner product of the reaction field factors f_ll'mm'(J,K)
!!    and the multipole moments M_lm(J) to obtain the reaction field
!!    coefficients R_l'm'(K).
!!
!!  USES
!   implicit none
!
!!  ARGUMENTS
!   real(dp), intent(in) :: point(:)
!   type(SolHarmBasisCenter), intent(in) :: expansion_centers(:)
!
!   real(dp), intent(out) :: distsqs(:)
!   real(dp), intent(out) :: relvecs(:,:)
!   integer, intent(out) :: sorted_indices(:)
!
!!  INPUTS
!!   o point -- coordinate of point of evaluation
!!   o expansion_centers -- expansion centers
!!  OUTPUT
!!   o distsqs -- squared distances of all centers
!!                to point of evaluation, unsorted
!!   o relvecs -- all relative vectors between expansion centers and
!!                point of evaluation, unsorted
!!   o sorted_indices -- list of indices, sorted for ascending order of
!!                       squared distance
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  HISTORY
!!    Development version, FHI-aims (2014).
!!  SOURCE
!
!   integer :: i_center
!
!   do i_center = lbound(expansion_centers,1), ubound(expansion_centers,1)
!      relvecs(:,i_center) = point - expansion_centers(i_center)%coord
!      distsqs(i_center) = dot_product(relvecs(:,i_center), relvecs(:,i_center))
!      sorted_indices(i_center) = i_center ! save index
!   enddo
!
!   call dquicksort_indexlist(size(expansion_centers), distsqs, sorted_indices)
!
!end subroutine sort_centers
!!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/center_convolution
!  NAME
!    center_convolution
!  SYNOPSIS

pure function center_convolution(center, numbasis) result(val)

!  PURPOSE
!    Calculate the inner product of coefficients and numerically evaluated
!    basis functions.
!
!  USES
   implicit none

!  ARGUMENTS
   type(SolHarmBasisCenter), intent(in) :: center
   real(dp), intent(in) :: numbasis(:)

   real(dp) :: val
!  INPUTS
!   o center -- has sparse expansion coefficients attached
!   o numbasis -- spherical harmonics basis evaluated at a certain point
!  OUTPUT
!   o val -- value of the potential
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: i_sp, ilm0

   val = 0.e0_dp
   ilm0 = SolHarmBasisCenter_lbound_lm(center) - 1
   do i_sp = lbound(center%coeff,1), ubound(center%coeff,1)
      val = val + center%coeff(i_sp)%val * &
                     numbasis(center%coeff(i_sp)%ind + ilm0)
   enddo ! i_sp

end function center_convolution
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/calculate_solvent_potential_at_single_point
!  NAME
!    calculate_solvent_potential_at_single_point
!  SYNOPSIS

subroutine calculate_potential_at_points(n_p, points, b, potential)

!  PURPOSE
!    Evaluate the reaction field at a given point at several expansion
!    centers and interpolate between the results.
!
!  USES
   implicit none
!TODO change to fixed dimensions for now in order to treat single points
!  ARGUMENTS
   integer, intent(in) :: n_p
   real(dp), intent(in) :: points(3,n_p)
   class(Basis), intent(in) :: b

   real(dp), intent(out) :: potential(n_p)

!  INPUTS
!   o points -- coordinate of points of evaluation
!   o b -- basis with expansion coefficients
!  OUTPUT
!   o potential -- interpolated value of the reaction field
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'calculate_potential_at_points'
   character(132) :: info_str
   integer :: max_lmax, max_lm, n_centers, i_point
   real(dp), allocatable :: cartesians(:,:), rlylm(:)

   select type (b)
   class is (SolHarmBasis)
      max_lmax = maxval(b%centers(:)%lmax)
      max_lm = (max_lmax+1)**2
      n_centers = size(b%centers)

      ! initialize cartesians
      call initialize_cartesian_ylm(max_lmax)

      call f_allocate(cartesians, 1,n_max_cartesian, 0,max_lmax, name="cartesians")
      call f_allocate(rlylm, max_lm, name="rlylm")

      ! loop over points
      select case (b%solharm_type)
      case (MPE_CONST%BASIS_IRR)
         do i_point = 1, size(points,2)
            call calculate_irregular_at_point( points(:,i_point), &
                  b, potential(i_point), &
                  cartesians, rlylm )
         enddo ! i_point

      case (MPE_CONST%BASIS_REG)
         do i_point = 1, size(points,2)
            call calculate_regular_at_point( points(:,i_point), &
                  b, potential(i_point), &
                  cartesians, rlylm )
         enddo ! i_point

      case default
         write(info_str,'(A,I2,A)') '*** ERROR: basis type no. ', &
                              b%solharm_type, ' is unknown'
         call localorb_info(info_str)
         call f_stop('Encountered unknown basis type', func)
      endselect

      call f_deallocate(rlylm, name="rlylm")
      call f_deallocate(cartesians, name="cartesians")

   class default
      call f_stop('Internal Error: Potential evaluation not implemented '//&
                     ' for this basis type', func)
   end select

   contains


   subroutine calculate_irregular_at_point( point, b, potential, &
                                 cartesians, rlylm )
      implicit none

      real(dp), intent(in) :: point(:)
      type(SolHarmBasis), intent(in) :: b
      real(dp), intent(out) :: potential
      real(dp) :: cartesians(:,:), rlylm(:)

      integer :: i_center, n_centers, i_l, i_m, i_lm
      real(dp) :: relvec(3), rsq_inv, r_to_m2lm1

      potential = 0.e0_dp

      n_centers = size(b%centers)

      do i_center = 1, n_centers

         ! get distance vector
         relvec = point - b%centers(i_center)%coord

         ! calculate potential
         call evaluate_onecenter_cartesians( relvec, &
                              b%centers(i_center)%lmax, cartesians )
         call tab_ylm_onecenter_cartesian( b%centers(i_center)%lmax, &
                              cartesians, rlylm )

         rsq_inv = 1.e0_dp/dot_product(relvec, relvec)
         r_to_m2lm1 = sqrt(rsq_inv)
         i_lm = 0
         do i_l = 0, b%centers(i_center)%lmax
            do i_m = -i_l, i_l
               i_lm = i_lm + 1
               rlylm(i_lm) = rlylm(i_lm) * r_to_m2lm1
            enddo ! i_m
            r_to_m2lm1 = r_to_m2lm1 * rsq_inv
         enddo ! i_l

         potential = potential + &
                     center_convolution( b%centers(i_center), rlylm )
      enddo ! i_center

   end subroutine calculate_irregular_at_point

   subroutine calculate_regular_at_point( point, b, potential, &
                                 cartesians, rlylm )
      implicit none

      real(dp), intent(in) :: point(:)
      type(SolHarmBasis), intent(in) :: b
      real(dp), intent(out) :: potential
      real(dp) :: cartesians(:,:), rlylm(:)

      integer :: i_center, n_centers
      real(dp) :: relvec(3)

      potential = 0.e0_dp

      n_centers = size(b%centers)

      do i_center = 1, n_centers

         ! get distance vector
         relvec = point - b%centers(i_center)%coord

         ! calculate potential
         call evaluate_onecenter_cartesians( relvec, &
                              b%centers(i_center)%lmax, cartesians )
         call tab_ylm_onecenter_cartesian( b%centers(i_center)%lmax, &
                              cartesians, rlylm )
         potential = potential + &
                     center_convolution( b%centers(i_center), rlylm )
      enddo ! i_center

   end subroutine calculate_regular_at_point

end subroutine calculate_potential_at_points
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/isqrt_newton
!  NAME
!    isqrt_newton
!  SYNOPSIS
elemental function isqrt_newton(n) result(x)
!  USES
   implicit none
!  PURPOSE
!    Compute the integer square root of n using Newton's method
!  ARGUMENTS
   integer, intent(in) :: n
   integer :: x
!  INPUTS
!    o n -- integer number
!  RETURNS
!    o x -- integer square root of n
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE
   integer :: y
   x = n
   y = (x + 1) / 2
   do while (y < x)
      x = y
      y = (x + n/x) / 2
   enddo
end function isqrt_newton
!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/map_lm_to_l
!  NAME
!    map_lm_to_l
!  SYNOPSIS
elemental function map_lm_to_l(i_lm) result(i_l)
!  USES
   implicit none
!  PURPOSE
!    Restore the l value from the combined l,m index
!  ARGUMENTS

   integer, intent(in) :: i_lm
   integer :: i_l

!  INPUTS
!    o i_lm -- combined l,m index
!  OUTPUT
!    o i_l -- l quantum number
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE
   i_l = isqrt_newton(i_lm)
end function map_lm_to_l
!******

end module mpe_reaction_field

