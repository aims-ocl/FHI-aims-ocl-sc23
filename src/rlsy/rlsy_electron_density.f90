module rlsy_electron_density
!!
!! Someplace to store the electron density in the irreducible representation.
!!
use rlsy_constants, only: r8,r16,rl_huge,rl_hugeint,rl_twopi
use rlsy_helpers, only: tochar,rl_chop,norm2
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully
use rlsy_integration_grid, only: rl_integration_grid
use rlsy_free_atom, only: rl_free_atom
use rlsy_extended_cluster, only: rl_extended_cluster
use rlsy_crystalstructure, only: rl_crystalstructure

implicit none
private
public :: rl_electron_density

!> container that holds density. Really not convinced this is a good idea at all.
type rl_electron_density
    !> how many spin channels
    integer :: n_spin=-rl_hugeint
    !> how many electrons should it be
    real(r8) :: n_electron=-rl_huge
    !> how many semilocal integration points
    integer :: n_semilocal_point=-rl_hugeint
    !> normal electron density (n_semilocal_point,n_spin)
    real(r8), dimension(:,:), allocatable :: irr_rho
    !> normal electron density gradient (3,n_semilocal_point,n_spin)
    real(r8), dimension(:,:,:), allocatable :: irr_grad_rho
    !> change in electron density (n_semilocal_point,n_spin)
    real(r8), dimension(:,:), allocatable :: irr_delta_rho

    !> irreducible free atom density (n_semilocal_point)
    real(r8), dimension(:), allocatable :: irr_free_rho
    !> irreducible free atom density gradient (3,n_semilocal_point)
    real(r8), dimension(:,:), allocatable :: irr_free_grad_rho

    contains
        !> initialize storage
        procedure :: generate
end type

contains

!> initialize storage for irreducible charge density
subroutine generate(density,grid,p,free,ec,n_spin,n_electron,mw)
    !> electron density
    class(rl_electron_density), intent(out) :: density
    !> integration grid
    type(rl_integration_grid), intent(inout) :: grid
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> free atom properties
    type(rl_free_atom), intent(in) :: free
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> number of spins
    integer, intent(in) :: n_spin
    !> number of electrons
    real(r8), intent(in) :: n_electron
    !> mpi helper
    type(rl_mpi_helper), intent(inout) :: mw
    !init: block
    integer :: ib
    !freeatom: block
    real(r8), dimension(3) :: v
    real(r8) :: total_charge,scalefactor
    real(r8) :: r
    !integer :: i,ib,ip,jp,ie,ispc,iatm
    integer :: i,ip,jp,ie,ispc,iatm

    ! set some basic things
    !init: block
        ! Number of spin channels
        density%n_spin = n_spin
        ! Number of electrons
        density%n_electron = n_electron

        ! Count number of semilocal points. These are the
        ! number of irreducible points in this chunk.
        density%n_semilocal_point=0
        do ib=1,grid%n_irr_batch
            density%n_semilocal_point=density%n_semilocal_point+grid%irr_batch(ib)%n_point
        enddo
        call grid%ml%allreduce('sum',density%n_semilocal_point)

        ! Make some space
        allocate(density%irr_free_rho(density%n_semilocal_point))
        allocate(density%irr_free_grad_rho(3,density%n_semilocal_point))
        allocate(density%irr_rho(density%n_semilocal_point,density%n_spin))
        allocate(density%irr_grad_rho(3,density%n_semilocal_point,density%n_spin))
        allocate(density%irr_delta_rho(density%n_semilocal_point,density%n_spin))
        density%irr_free_rho=0.0_r8
        density%irr_free_grad_rho=0.0_r8
        density%irr_rho=0.0_r8
        density%irr_grad_rho=0.0_r8
        density%irr_delta_rho=0.0_r8
    !end block init

    ! Calculate the free atom density just to make sure that I can, not for any particular purpose.
    !freeatom: block

        total_charge=0.0_r8

        do ib=1,grid%n_irr_batch
        do ip=1,grid%irr_batch(ib)%n_point
            ! Index in rho
            jp=grid%irr_batch(ib)%semilocal_irr_offset+ip
            ! Go through extended atoms
            do i=1,grid%irr_batch(ib)%n_relevant_ext_atom
                ! Index to extended cluster
                ie=grid%irr_batch(ib)%relevant_ext_atom(i)
                iatm=ec%index_unit_cell(ie)
                ispc=p%species(iatm)
                ! Distance to integration points
                v=ec%cartesian_coordinate(:,ie) - grid%irr_batch(ib)%folded_coordinate(:,ip)
                r=norm2( v )
                ! Accumulate density
                density%irr_free_rho(jp)=density%irr_free_rho(jp)+free%species(ispc)%rho%evaluate(r)
                ! And density gradient
                if ( r .gt. 1E-10_r8 ) then
                    ! Note negative sign here to get the Cartesian gradient from the scalar derivative.
                    v=-free%species(ispc)%rho_deriv%evaluate(r)*v/r
                else
                    v=0.0_r8
                endif
                density%irr_free_grad_rho(:,jp)=density%irr_free_grad_rho(:,jp)+v
            enddo
            ! Accumulate total charge
            total_charge=total_charge+density%irr_free_rho(jp)*grid%irr_batch(ib)%integration_weight(ip)*grid%irr_batch(ib)%partition_function(ip)
        enddo
        enddo
        ! Sync density
        call grid%ml%allreduce('sum',density%irr_free_rho)
        call grid%ml%allreduce('sum',density%irr_free_grad_rho)
        ! Sync total charge
        call mw%allreduce('sum',total_charge)
        ! For some reason there is a factor 4*pi between the free atom density and the normal one.
        ! I let the scaling factor take care of that. Will have to check the gradient later.
        !total_charge=total_charge/(2.0_r8*rl_twopi)
        scalefactor=n_electron/total_charge
        density%irr_free_rho=density%irr_free_rho*scalefactor
        density%irr_free_grad_rho=density%irr_free_grad_rho*scalefactor
    !end block freeatom

end subroutine

end module
