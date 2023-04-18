module rlsy_realspace_integration_helper
!!
!! Thing that holds temporary storage during integrations. Made this its own thing
!! to avoid boilerplate copy-paste all over the place.
!!
use rlsy_constants, only: r8,rl_iou,rl_exitcode_param,rl_exitcode_memory,rl_hugeint,rl_pi
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_helpers, only: tochar,rl_mom_real,norm2
use rlsy_sorting, only: rl_qsort
use rlsy_integration_grid, only: rl_integration_grid,rl_integration_batch
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_extended_cluster, only: rl_extended_cluster,rl_extended_cluster_hashedspace
use rlsy_realspace_matrix, only: rl_realspace_matrix,rl_realspace_matrix_notdistributed,&
    rl_realspace_matrix_mediumdistributed,rl_realspace_matrix_fulldistributed
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc,&
    rl_kspace_eigenproblem_prob,rl_kspace_eigenproblem_prob_real,&
    rl_kspace_eigenproblem_prob_complex,rl_kspace_eigenproblem_kpoint_real,&
    rl_kspace_eigenproblem_kpoint_complex
use rlsy_verletlist, only: rl_verletbox
use rlsy_electron_density, only: rl_electron_density
use rlsy_timer, only: rl_timer_density
use rlsy_pair_locator, only: rl_pair_locator

implicit none
private

public :: rl_integration_helper
public :: rl_integration_helper_multipoleexpansion

!> Generic integration helper with only the bare minimum. Holds temporary storage for one batch och integration points.
type, abstract :: rl_integration_helper
    !> some generic settings:
    integer :: derivative_level=-rl_hugeint

    !> how many active atoms are there (active means it can interact with some point in the batch)
    integer :: n_active_atom=-rl_hugeint
    !> list of active atoms. Indices points to atoms in the extended cluster
    integer, dimension(:), allocatable :: active_atom

    !> how many integration points are we working with
    integer :: n_integration_point=-rl_hugeint
    !> tabulated coordinates (n_active_atom,n_integration_point)
    !> by coordinate I mean the vectors that go from active atoms to integration points.
    !> separated into a direction vector (x,y,z), and a radial distance (r), norm2([x,y,z])=1
    real(r8), dimension(:,:), allocatable :: x  !< x/norm(r)
    real(r8), dimension(:,:), allocatable :: y  !< y/norm(r)
    real(r8), dimension(:,:), allocatable :: z  !< z/norm(r)
    real(r8), dimension(:,:), allocatable :: r  !< norm(r)
    real(r8), dimension(:,:), allocatable :: ir !< 1/norm(r)
    !> basis function index things
    integer :: n_nonzero_basis=-rl_hugeint
    integer, dimension(:), allocatable :: active_atom_offset
    integer, dimension(:), allocatable :: active_atom_nbasis
    integer, dimension(:,:), allocatable :: basis_ind

    !> helper arrays that hold intermediate basis function values
    real(r8), dimension(:,:), allocatable :: f          !< radial function
    real(r8), dimension(:,:), allocatable :: df         !< radial derivative
    real(r8), dimension(:,:), allocatable :: ddf        !< radial second derivative
    real(r8), dimension(:,:), allocatable :: Ylm        !< Ylm spherical harmonic
    real(r8), dimension(:,:), allocatable :: dYlmdx     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: dYlmdy     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: dYlmdz     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: ddYlmdxdx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdxdy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdxdz  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydz  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdz  !< Ylm Hessian

    ! !> flat arrays that hold basis function values and derivatives
    ! real(r8), dimension(:), allocatable :: basis_V          !< basis fn
    ! real(r8), dimension(:), allocatable :: basis_dVdx       !< basis gradient
    ! real(r8), dimension(:), allocatable :: basis_dVdy       !< basis gradient
    ! real(r8), dimension(:), allocatable :: basis_dVdz       !< basis gradient
    ! real(r8), dimension(:), allocatable :: basis_ddVdxdx    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdydx    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdzdx    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdxdy    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdydy    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdzdy    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdxdz    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdydz    !< basis Hessian
    ! real(r8), dimension(:), allocatable :: basis_ddVdzdz    !< basis Hessian

    ! real(r8), dimension(:,:,:), allocatable :: dm_D
    ! real(r8), dimension(:,:), allocatable :: dm_DtimesPhi
    ! real(r8), dimension(:,:), allocatable :: dm_Phi
    ! real(r8), dimension(:,:), allocatable :: dm_Phi_gx
    ! real(r8), dimension(:,:), allocatable :: dm_Phi_gy
    ! real(r8), dimension(:,:), allocatable :: dm_Phi_gz
    ! real(r8), dimension(:,:), allocatable :: dm_PhiT
    ! real(r8), dimension(:,:), allocatable :: dm_PhiT_gx
    ! real(r8), dimension(:,:), allocatable :: dm_PhiT_gy
    ! real(r8), dimension(:,:), allocatable :: dm_PhiT_gz

    ! !> density matrix with everything except the relevant entries zeroed out
    ! real(r8), dimension(:,:), allocatable :: pruned_density_matrix

    !contains
        ! !> initialize temporary storage
        ! procedure :: setup
        ! !> reset the work arrays in the active set to nothing.
        ! procedure :: reset
        ! !> update list of active basis functions
        ! procedure :: update_active_basis_functions
        ! !> tabulate coordinates
        ! procedure :: tabulate_coordinates
        ! !> tabulate basis functions
        ! procedure :: tabulate_basis_functions
        ! !> fetch pruned density matrix
        ! procedure :: fetch_pruned_density_matrix
        ! !> evaluate density
        ! procedure :: evaluate_density
end type

!> helper for density-like integrations
type, extends(rl_integration_helper) :: rl_integration_helper_density
    !> sensible to know the number of spin channels
    !integer :: n_spin=-rl_hugeint
    ! !> density
    ! real(r8), dimension(:,:), allocatable :: rho
    ! real(r8), dimension(:,:,:), allocatable :: grad_rho
    ! real(r8), dimension(:,:,:), allocatable :: hess_rho
end type

!> helper for overlap/hamiltonian-like integrations
type, extends(rl_integration_helper) :: rl_integration_helper_hamiltonian
    !> sensible to know the number of spin channels
    !integer :: n_spin=-rl_hugeint
    ! !> density
    ! real(r8), dimension(:,:), allocatable :: rho
    ! real(r8), dimension(:,:,:), allocatable :: grad_rho
    ! real(r8), dimension(:,:,:), allocatable :: hess_rho
end type

!> helper for determining multipole expansion coefficients.
type :: rl_integration_helper_multipoleexpansion
    !> Largest angular index (=l^2+l+1+m)
    integer :: n_iang=-rl_hugeint
    !> Current number of integration points
    integer :: n_integration_point=-rl_hugeint
    !> Number of radial shells
    integer :: n_radial_shell=-rl_hugeint

    !> Flattened arrays of integration points
    real(r8), dimension(:), allocatable :: x  !< x/norm(r)
    real(r8), dimension(:), allocatable :: y  !< y/norm(r)
    real(r8), dimension(:), allocatable :: z  !< z/norm(r)
    real(r8), dimension(:), allocatable :: r  !< norm(r)
    real(r8), dimension(:), allocatable :: ir !< 1/norm(r)
    !> integration weight
    real(r8), dimension(:), allocatable :: w
    !> delta-density, rho-rho_free
    real(r8), dimension(:), allocatable :: delta_rho
    !> delta-densitygradient, rho-rho_free, radial direction
    real(r8), dimension(:), allocatable :: delta_rho_grad
    !> spherical harmonics
    real(r8), dimension(:,:), allocatable :: Ylm
    !> buffer to store expansion (iang,irad,iatom)
    real(r8), dimension(:,:,:), allocatable :: coeffbuf
    !> buffer to store radial values
    real(r8), dimension(:,:), allocatable :: rvalbuf

    !> atom index, per point
    integer, dimension(:), allocatable :: ind_atm
    !> species index, per point
    !integer, dimension(:), allocatable :: ind_spc
    !> radial index, per point
    integer, dimension(:), allocatable :: ind_rad
    !> angular index, per point
    !integer, dimension(:), allocatable :: ind_ang
end type

contains

! !> tabulate the coordinates of all atoms with respect to all points
! subroutine tabulate_coordinates(active,ec,batch)
!     !> subset of basis functions and atoms
!     class(rl_active_set), intent(inout) :: active
!     !> extended cluster of atoms
!     type(rl_extended_cluster), intent(in) :: ec
!     !> batch information
!     class(rl_integration_batch), intent(in) :: batch
!
!     real(r8), dimension(3) :: v,w
!     real(r8) :: f0,f1,tol
!     integer :: np,na,i,j,k
!
!     ! Some shorthand
!     np=batch%n_point
!     na=active%n_active_atom
!     tol=1E-10_r8
!
!     ! calculate coordinates
!     do i=1,na
!         k=active%active_atom(i)
!         w=ec%cartesian_coordinate(:,k)
!         do j=1,np
!             v=batch%folded_coordinate(:,j)-w
!             f0=norm2(v)
!             if ( f0 .gt. tol ) then
!                 f1=1.0_r8/f0
!                 v=v*f1
!                 active%x(j,i)=v(1)
!                 active%y(j,i)=v(2)
!                 active%z(j,i)=v(3)
!                 active%r(j,i)=f0
!                 active%ir(j,i)=f1
!             else
!                 active%x(j,i)=0.0_r8
!                 active%y(j,i)=0.0_r8
!                 active%z(j,i)=0.0_r8
!                 active%r(j,i)=f0
!                 active%ir(j,i)=0.0_r8
!             endif
!         enddo
!     enddo
! end subroutine

! !> convert the list of active atoms to a list of active basis functions
! subroutine update_active_basis_functions(active,ec,p,basis)
!     !> subset of basis functions and atoms
!     class(rl_active_set), intent(inout) :: active
!     !> extended cluster of atoms
!     type(rl_extended_cluster), intent(in) :: ec
!     !> structure
!     type(rl_crystalstructure), intent(in) :: p
!     !> basis set information
!     type(rl_lcao_basis_set), intent(in) :: basis
!
!
!     integer :: i,j,l
!     integer :: is,iu,ie,nb
!
!     ! This could probably be made faster by first generating a possible pool for
!     ! all the batches or something. Don't think this is a bottleneck though.
!     active%basis_ind=0
!     l=0
!     do i=1,active%n_active_atom
!         ie=active%active_atom(i)
!         iu=ec%index_unit_cell(ie)
!         is=p%species(iu)
!         nb=basis%species(is)%n_basis
!         active%active_atom_offset(i)=l
!         active%active_atom_nbasis(i)=nb
!         do j=1,nb
!             l=l+1
!             active%basis_ind(:,l)=[basis%species(is)%basis_fn(j),basis%species(is)%basis_ang(j)]
!         enddo
!     enddo
!     active%n_nonzero_basis=l
!
!     ! Equipped with this information, I can make some temporary space:
!     select case(active%derivative_level)
!     case(0)
!         allocate(active%dm_D(active%n_nonzero_basis,active%n_nonzero_basis,active%n_spin))
!         allocate(active%dm_DtimesPhi(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_Phi(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_PhiT(active%n_integration_point,active%n_nonzero_basis))
!     case(1)
!         allocate(active%dm_D(active%n_nonzero_basis,active%n_nonzero_basis,active%n_spin))
!         allocate(active%dm_DtimesPhi(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_Phi(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_Phi_gx(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_Phi_gy(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_Phi_gz(active%n_nonzero_basis,active%n_integration_point))
!         allocate(active%dm_PhiT(active%n_integration_point,active%n_nonzero_basis))
!         allocate(active%dm_PhiT_gx(active%n_integration_point,active%n_nonzero_basis))
!         allocate(active%dm_PhiT_gy(active%n_integration_point,active%n_nonzero_basis))
!         allocate(active%dm_PhiT_gz(active%n_integration_point,active%n_nonzero_basis))
!     end select
! end subroutine

! !> tabulate the values of the wave functions over the batch
! subroutine tabulate_basis_functions(active,basis)
!     !> subset of basis functions and atoms
!     class(rl_active_set), intent(inout) :: active
!     !> basis set information
!     type(rl_lcao_basis_set), intent(in) :: basis
!
!     integer :: i,j,i1,i2,ir,ia,np
!
!     np=active%n_integration_point
!     ! Do it slightly different depending on how many derivatives we need.
!     select case(active%derivative_level)
!     case(0)
!         ! No derivatives, just function values
!         do i=1,active%n_active_atom
!             ! Evaluate radial functions
!             do ir=1,basis%n_radial
!                 call basis%radial(ir)%val( active%r(1:np,i),active%f(1:np,ir) )
!             enddo
!             ! Evaluate spherical harmonics
!             do ia=1,basis%n_angular
!                 call basis%Ylm( active%x(1:np,i),active%y(1:np,i),active%z(1:np,i),ia,active%Ylm(1:np,ia) )
!             enddo
!             ! Combine basis functions to the full ones
!             do j=active%active_atom_offset(i)+1,active%active_atom_offset(i)+active%active_atom_nbasis(i)
!                 ir = active%basis_ind(1,j) ! radial index of basis fn j
!                 ia = active%basis_ind(2,j) ! angular index of basis fn j
!                 ! where to store values in flat array
!                 i1=(j-1)*np+1
!                 i2=j*np
!                 !active%basis_v(i1:i2)=active%f(1:np,ir) * active%Ylm(1:np,ia)
!                 active%dm_PhiT(:,j)=active%f(1:np,ir) * active%Ylm(1:np,ia)
!             enddo
!         enddo
!         active%dm_Phi=transpose(active%dm_PhiT)
!     case(1)
!         ! Get values + gradient
!         do i=1,active%n_active_atom
!
!             ! Evaluate radial functions + first derivative
!             do ir=1,basis%n_radial
!                 call basis%radial(ir)%val_der( active%r(1:np,i),active%f(1:np,ir),active%df(1:np,ir) )
!             enddo
!             ! Evaluate spherical harmonics + gradient
!             do ia=1,basis%n_angular
!                 call basis%Ylm_grad( active%x(1:np,i),active%y(1:np,i),active%z(1:np,i),active%ir(1:np,i),ia,&
!                 active%Ylm(1:np,ia),active%dYlmdx(1:np,ia),active%dYlmdy(1:np,ia),active%dYlmdz(1:np,ia) )
!             enddo
!
!             ! Loop over the basis functions on this atom.
!             do j=active%active_atom_offset(i)+1,active%active_atom_offset(i)+active%active_atom_nbasis(i)
!                 ir = active%basis_ind(1,j) ! radial index of basis fn j
!                 ia = active%basis_ind(2,j) ! angular index of basis fn j
!                 ! where to store values in flat array
!                 i1=(j-1)*np+1
!                 i2=j*np
!                 ! Store function values
!                 active%dm_PhiT(:,j)=active%f(1:np,ir) * active%Ylm(1:np,ia)
!                 ! Store gradient
!                 active%dm_PhiT_gx(:,j)= active%x(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdx(1:np,ia)
!                 active%dm_PhiT_gy(:,j)= active%y(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdy(1:np,ia)
!                 active%dm_PhiT_gz(:,j)= active%z(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdz(1:np,ia)
!                 ! active%dm_PhiT_gx(:,j)= active%dYlmdx(1:np,ia)
!                 ! active%dm_PhiT_gy(:,j)= active%dYlmdy(1:np,ia)
!                 ! active%dm_PhiT_gz(:,j)= active%dYlmdz(1:np,ia)
!
!             enddo
!         enddo
!         active%dm_Phi=transpose(active%dm_PhiT)
!         active%dm_Phi_gx=transpose(active%dm_PhiT_gx)
!         active%dm_Phi_gy=transpose(active%dm_PhiT_gy)
!         active%dm_Phi_gz=transpose(active%dm_PhiT_gz)
!     case(2)
!         ! Get values + gradient + hessian
!     end select
! end subroutine

! !> reset the active set for work on a new batch
! subroutine reset(active,ec,basis,grid,batch,mem)
!     !> subset of basis functions and atoms
!     class(rl_active_set), intent(inout) :: active
!     !> extended cluster of atoms
!     type(rl_extended_cluster), intent(in) :: ec
!     !> basis set information
!     type(rl_lcao_basis_set), intent(in) :: basis
!     !> grid information
!     type(rl_integration_grid), intent(in) :: grid
!     !> batch information
!     class(rl_integration_batch), intent(in) :: batch
!     !> memory tracker
!     type(rl_memtracker), intent(in) :: mem
!
!     ! Zero all basis functions
!     active%n_nonzero_basis    = 0
!     active%active_atom_offset = 0
!     active%active_atom_nbasis = 0
!     active%basis_ind          = 0
!
!     ! Make a note of the number of integration points
!     active%n_integration_point=batch%n_point
!
!     ! set the active atoms to the ones determined by the batch?
!     active%n_active_atom=batch%n_relevant_ext_atom
!     active%active_atom = 0
!     active%active_atom(1:active%n_active_atom)=batch%relevant_ext_atom(1:active%n_active_atom)
!
!     ! space for coordinates
!     active%x  = 0.0_r8
!     active%y  = 0.0_r8
!     active%z  = 0.0_r8
!     active%r  = 0.0_r8
!     active%ir = 0.0_r8
!     select case(active%derivative_level)
!     case(0)
!         active%rho     = 0.0_r8
!         ! space for basis function values
!         active%basis_V = 0.0_r8
!         active%f       = 0.0_r8
!         active%Ylm     = 0.0_r8
!     case(1)
!         active%rho     = 0.0_r8
!         active%grad_rho= 0.0_r8
!         ! space for basis function values
!         active%basis_V    = 0.0_r8
!         active%basis_dVdx = 0.0_r8
!         active%basis_dVdy = 0.0_r8
!         active%basis_dVdz = 0.0_r8
!         ! Intermediate things for derivatives
!         active%f      = 0.0_r8
!         active%df     = 0.0_r8
!         active%Ylm    = 0.0_r8
!         active%dYlmdx = 0.0_r8
!         active%dYlmdy = 0.0_r8
!         active%dYlmdz = 0.0_r8
!     case(2)
!         ! space for basis function values
!         active%basis_V       = 0.0_r8
!         active%basis_dVdx    = 0.0_r8
!         active%basis_dVdy    = 0.0_r8
!         active%basis_dVdz    = 0.0_r8
!         active%basis_ddVdxdx = 0.0_r8
!         active%basis_ddVdydx = 0.0_r8
!         active%basis_ddVdzdx = 0.0_r8
!         active%basis_ddVdxdy = 0.0_r8
!         active%basis_ddVdydy = 0.0_r8
!         active%basis_ddVdzdy = 0.0_r8
!         active%basis_ddVdxdz = 0.0_r8
!         active%basis_ddVdydz = 0.0_r8
!         active%basis_ddVdzdz = 0.0_r8
!         ! Intermediate things for derivatives
!         active%f         = 0.0_r8
!         active%df        = 0.0_r8
!         active%ddf       = 0.0_r8
!         active%Ylm       = 0.0_r8
!         active%dYlmdx    = 0.0_r8
!         active%dYlmdy    = 0.0_r8
!         active%dYlmdz    = 0.0_r8
!         active%ddYlmdxdx = 0.0_r8
!         active%ddYlmdydx = 0.0_r8
!         active%ddYlmdzdx = 0.0_r8
!         active%ddYlmdxdy = 0.0_r8
!         active%ddYlmdydy = 0.0_r8
!         active%ddYlmdzdy = 0.0_r8
!         active%ddYlmdxdz = 0.0_r8
!         active%ddYlmdydz = 0.0_r8
!         active%ddYlmdzdz = 0.0_r8
!     end select
! end subroutine
!
! !> make some temporary space
! subroutine setup(active,ec,basis,grid,mem,n_spin,derivative_level)
!     !> subset of basis functions and atoms
!     class(rl_active_set), intent(inout) :: active
!     !> extended cluster of atoms
!     type(rl_extended_cluster), intent(in) :: ec
!     !> basis set information
!     type(rl_lcao_basis_set), intent(in) :: basis
!     !> grid information
!     type(rl_integration_grid), intent(in) :: grid
!     !> memory tracker
!     type(rl_memtracker), intent(in) :: mem
!     !> number of spin channels
!     integer, intent(in) :: n_spin
!     !> how many derivatives are to be taken
!     integer, intent(in) :: derivative_level
!
!     ! Set some constants
!     active%derivative_level=derivative_level
!     active%n_spin = n_spin
!
!     ! Space for tracking relevant atoms and basis functions
!     allocate(active%active_atom(        grid%max_n_atom_per_batch) )
!     allocate(active%active_atom_offset( grid%max_n_atom_per_batch) )
!     allocate(active%active_atom_nbasis( grid%max_n_atom_per_batch) )
!     allocate(active%basis_ind(        2,grid%max_n_basis_per_batch))
!
!     ! Space for coordinates
!     allocate(active%x (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
!     allocate(active%y (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
!     allocate(active%z (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
!     allocate(active%r (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
!     allocate(active%ir(grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
!
!
!     ! space for pruned density matrix @TODO don't forget nspin here, eventually.
!     allocate(active%pruned_density_matrix(grid%max_n_basis_per_batch**2,active%n_spin))
!
!     select case(active%derivative_level)
!     case(0)
!         ! Space for density
!         allocate(active%rho(grid%max_n_point_per_batch,active%n_spin))
!         ! space for basis function values
!         allocate(active%basis_V      (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         ! Intermediate things for derivatives
!         allocate(active%f      ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%Ylm    ( grid%max_n_point_per_batch,basis%n_angular ))
!     case(1)
!         ! Space for density
!         allocate(active%rho(grid%max_n_point_per_batch,active%n_spin))
!         ! Space for gradient
!         allocate(active%grad_rho(3,grid%max_n_point_per_batch,active%n_spin))
!         ! space for basis function values
!         allocate(active%basis_V      (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdx   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdy   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdz   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         ! Intermediate things for derivatives
!         allocate(active%f      ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%df     ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%Ylm    ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdx ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdy ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdz ( grid%max_n_point_per_batch,basis%n_angular ))
!     case(2)
!         ! space for basis function values
!         allocate(active%basis_V       (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdx    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdy    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_dVdz    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdxdx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdydx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdzdx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdxdy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdydy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdzdy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdxdz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdydz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         allocate(active%basis_ddVdzdz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
!         ! Intermediate things for derivatives
!         allocate(active%f         ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%df        ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%ddf       ( grid%max_n_point_per_batch,basis%n_radial  ))
!         allocate(active%Ylm       ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdx    ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdy    ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%dYlmdz    ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdxdx ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdydx ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdzdx ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdxdy ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdydy ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdzdy ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdxdz ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdydz ( grid%max_n_point_per_batch,basis%n_angular ))
!         allocate(active%ddYlmdzdz ( grid%max_n_point_per_batch,basis%n_angular ))
!     end select
! end subroutine

end module
