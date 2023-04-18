module rlsy_hartree_potential
!!
!! Handles everyting electrostatic
!!
use rlsy_constants, only: r8,rl_iou,rl_huge,rl_tiny,rl_hugeint,rl_exitcode_memory,rl_twopi
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
!use rlsy_helpers, only: tochar,rl_mom_real
!use rlsy_sorting, only: rl_qsort
use rlsy_helpers, only: norm2
use rlsy_integration_grid, only: rl_integration_grid,rl_integration_batch!,rl_integration_batch_irr,rl_integration_batch_full
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_free_atom, only: rl_free_atom
use rlsy_extended_cluster, only: rl_extended_cluster,rl_extended_cluster_hashedspace
use rlsy_verletlist, only: rl_verletbox
use rlsy_electron_density, only: rl_electron_density
use rlsy_pair_locator, only: rl_pair_locator
use rlsy_realspace_integration_helper, only: rl_integration_helper_multipoleexpansion
use rlsy_timer, only: rl_timer_multipole

implicit none
private

public :: rl_multipole_expansion
public :: rl_calculate_multipole_expansion


!> spline in multipole expansion
type rl_hartree_multipole_spline
    !> smallest r (below this I switch to value that tapers to zero correctly)
    real(r8), private :: r_min=-rl_huge
    !> largest r (above this we use the far-field potential)
    real(r8), private :: r_max=-rl_huge
    !> number of knots
    integer, private :: n_knot=-rl_hugeint
    !> grid scaling
    real(r8), private :: rscale=-rl_hugeint
    !> l angular index for this spline
    integer, private :: l=-rl_hugeint
    !> fair-field multipole component
    real(r8), private :: farfield=-rl_huge
    !> pre-computed parameters to convert coordinate to indices
    real(r8), private :: invnsq=-rl_huge !< 1/(n_grid+1)^2
    real(r8), private :: np1=-rl_huge    !< n_grid + 1
    real(r8), private :: invC=-rl_huge   !< 1/rscale
    !> rho spline coefficients for very small r:
    real(r8), private :: sa2=-rl_huge
    real(r8), private :: sa3=-rl_huge
    real(r8), private :: sa4=-rl_huge
    !> density spline coefficients (spline is for rho(r)*r^2)
    real(r8), dimension(:,:), allocatable, private :: coeff_density
    !> potential 'spline' coefficients. Not really a spline.
    real(r8), dimension(:,:), allocatable, private :: coeff_potential
    contains
        procedure :: rho_val=>evaluate_rho_spline
        procedure :: rho_val_deriv=>evaluate_rho_spline_deriv
        !procedure :: potential_val=>evaluate_potential_function
        !procedure :: potential_val_deriv=>evaluate_potential_function_deriv
end type

!> atom in multipole expansion
type rl_multipole_expansion_atom_sym
    !> max l on this atom
    integer :: l_hartree
    !> number of multipole components
end type
!> atom in multipole expansion
type rl_multipole_expansion_atom_nonsym
    !> max l on this atom
    integer :: l_hartree
    !> number of multipole components
end type

!> multipole expansion of delta-density
type rl_multipole_expansion
    !> list of splines, may or may not be distributed.

    contains
        !> setup the handle
        procedure :: generate=>setup_multipole_expansion
end type

!@TODO Perhaps make it distributed in various ways.

contains

!> generate and set up the multipole expansion.
subroutine setup_multipole_expansion(mlt,l_hartree,p,sym,grid,mw,mem,verbosity)
    !> handle for mulitpole expansion
    class(rl_multipole_expansion), intent(out) :: mlt
    !> max l per species (from AIMS)
    integer, dimension(:), intent(in) :: l_hartree
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: t0,t1,timer

    !symsetup: block
    real(r8), dimension(:,:), allocatable :: coeff
    real(r8), dimension(:,:), allocatable :: mI,mR
    integer :: irratom,iatm,ispc,lmax
    integer :: nang
    integer :: l

    ! Set some basic things
    !init: block
        ! Start timers
        timer=mpi_wtime()
        t0=timer
        t1=timer

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'SETTING UP MULTIPOLE EXPANSION'
        endif

        ! First make note of the largest l in the expansion

    !end block init

    ! Depending on wether we have symmetries or not, we will do it slightly differently.
    select case(sym%n_operation)
    case(1)
    ! Only trivial symmetry, use all components in expansion.
    !nosymsetup: block
        !
    !end block nosymsetup
    case default
    ! More than one symmetry, think about symmetry.
    !symsetup: block

        ! For each irreducible atom, determine
        do irratom=1,sym%n_irreducible_atom
            iatm=sym%irr_to_all(irratom)
            ispc=p%species(iatm)
            lmax=l_hartree(ispc)
            ! So, now we know the largest l, which atom and which species we are looking at.
            ! Start creating the irreducible representation thing. First, how many Ylm can
            ! there be?
            nang=(2*lmax+1)**2
            allocate(coeff(nang,nang))
            ! For l=0, we always have a component.
            coeff=0.0_r8
            coeff(1,1)=1.0_r8
            ! Now check the other l's.
            do l=0,lmax
            enddo
        enddo

    !end block symsetup
    end select

write(*,*) 'stopstop setup multipole'
call mw%destroy()
stop

end subroutine

!>
subroutine rl_calculate_multipole_expansion(density,grid,p,sym,basis,free,ec,mw,mem,verbosity,tmr)
    !> density
    type(rl_electron_density), intent(in) :: density
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> free atom quantities
    type(rl_free_atom), intent(in) :: free
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timer
    type(rl_timer_multipole), intent(inout) :: tmr

    ! helper tings needed for te integration
    type(rl_integration_helper_multipoleexpansion) :: mih
    integer, parameter :: lmax=4
    real(r8), dimension(8) :: mellantid
    real(r8) :: timer,t0,t1
    !init: block
    integer :: a1,nb1,nb2,ipair
    !expandinmultipoles: block
    integer :: ibatch

    !init: block

        ! Start timers
        timer=mpi_wtime()
        t0=timer
        t1=timer
        mellantid=0.0_r8

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CALCULATING HARTREE POTENTIAL'
        endif

        ! Several ways to distribute this is if it becomes too memory
        ! greedy. Don't do all atoms at once, only do a few l,m per
        ! iteration are easy ways to make it completely distributed.
        ! later problem. Now I do everything at once. First create
        ! a little buffer space:
        call setup_integration_helper_multipoleexpansion(mih,grid,p,sym,lmax,mw,mem)
        t1=mpi_wtime(); mellantid(1)=mellantid(1)+t1-t0; t0=t1

        if ( verbosity .gt. 0 ) then
            write(*,*) '... initialized integration helper'
            write(*,*) 'pot zero',free%species(1)%potential_at_nucleus
        endif
    !end block init

    ! Expand density in atom-centered multipoles.
    !expandinmultipoles: block

        ! I'm not even convinced the batches help here, if this turns
        ! out to be slow I can flatten and vectorize everything. But
        ! I don't think this is a bottleneck right now.
        do ibatch=1,grid%n_irr_batch
            ! Reset, that means fetch new coordinates and delta-rho, and destroy data from previous iteration.
            call reset_integration_helper_multipoleexpansion(mih,grid,density,p,sym,grid%irr_batch(ibatch))
            t1=mpi_wtime(); mellantid(2)=mellantid(2)+t1-t0; t0=t1
            ! Start accumulating for the expansion coefficients.
            call accumulate_expansion(mih,basis)
            t1=mpi_wtime(); mellantid(3)=mellantid(3)+t1-t0; t0=t1
        enddo
        ! Sync the coefficient buffer across all ranks
        call mw%allreduce('sum',mih%coeffbuf)
        call mw%allreduce('max',mih%rvalbuf)

        ! Turn the coefficients into splines. Hmmm.
        call create_potential_spline(mih%coeffbuf,mih%rvalbuf,grid,free,p,mw,mem)
    !end block expandinmultipoles

    ! Cleanup the helper that determines multipole moments.
    call destroy_integration_helper_multipoleexpansion(mih,mem)

    if ( mw%talk ) then
        write(*,*) 'timings multipole expansion'
        write(*,*) '        init: ',mellantid(1)
        write(*,*) '       reset: ',mellantid(2)
        write(*,*) '        bfun: ',mellantid(3)
    endif

    ! And finally, sort out the timings
    mellantid(8)=mpi_wtime()-timer                ! total timer
    mellantid(7)=mellantid(8)-sum(mellantid(1:6)) ! time idle
!    tmr%init        = tmr%init +        mellantid(1)
!    tmr%idle        = tmr%idle +        mellantid(9)
!    tmr%total       = tmr%total +       mellantid(10)
!    tmr%fetchdm     = tmr%fetchdm +     mellantid(2)
!    tmr%reset       = tmr%reset +       mellantid(3)
!    tmr%tab_coord   = tmr%tab_coord +   mellantid(4)
!    tmr%prune_basis = tmr%prune_basis + mellantid(5)
!    tmr%tab_basis   = tmr%tab_basis +   mellantid(6)
!    tmr%prune_dm    = tmr%prune_dm +    mellantid(7)
!    tmr%matrixop    = tmr%matrixop +    mellantid(8)

    ! Check for stupidity
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

call mw%destroy()
stop

end subroutine

!> Start calculating actual expansion coefficients
subroutine accumulate_expansion(mih,basis)
    !> integration helper
    type(rl_integration_helper_multipoleexpansion), intent(inout) :: mih
    !> electron density
    type(rl_lcao_basis_set), intent(in) :: basis

    real(r8) :: f0,f1
    integer :: iang,ip,np,iatm,irad

    ! Evaluate spherical harmonics on the integration points
    np=mih%n_integration_point
    do iang=1,mih%n_iang
        call basis%Ylm( mih%x(1:np),mih%y(1:np),mih%z(1:np),iang,mih%Ylm(1:np,iang) )
        ! Multiply in weights and deltarho
        !mih%Ylm(1:np,iang)=mih%Ylm(1:np,iang)*mih%w(1:np)*mih%delta_rho(1:np)
        do ip=1,np
            iatm=mih%ind_atm(ip)
            irad=mih%ind_rad(ip)
            f0=mih%Ylm(ip,iang)*mih%w(ip)*mih%delta_rho(ip)
            f1=mih%Ylm(ip,iang)*mih%w(ip)*mih%delta_rho_grad(ip)
            mih%coeffbuf(iang,irad,iatm)=mih%coeffbuf(iang,irad,iatm)+f0
        enddo
    enddo

    ! Accumulate in the right place
    ! do ip=1,np
    !     iatm=mih%ind_atm(ip)
    !     irad=mih%ind_rad(ip)
    !     do iang=1,mih%n_iang
    !         mih%coeffbuf(iang,irad,iatm)=mih%coeffbuf(iang,irad,iatm)+mih%Ylm(ip,iang)
    !     enddo
    ! enddo
end subroutine

!> reset the integration helper
subroutine reset_integration_helper_multipoleexpansion(mih,grid,density,p,sym,batch)
    !> integration helper
    type(rl_integration_helper_multipoleexpansion), intent(inout) :: mih
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> electron density
    type(rl_electron_density), intent(in) :: density
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> batch
    type(rl_integration_batch), intent(in) :: batch

    real(r8), parameter :: weightprefactor=sqrt(2*rl_twopi) ! Normalization thing, confusing.
    real(r8), dimension(3) :: v0,v1
    real(r8) :: r,ir,w,dlrho,dlrhoderiv
    integer :: i,j,k,l,iop,ispc,irad,iang,jatm,jrad

    ! Set coordinates to nothing
    mih%x=0.0_r8
    mih%y=0.0_r8
    mih%z=0.0_r8
    mih%r=0.0_r8
    mih%w=0.0_r8
    mih%ind_atm=0
    mih%ind_rad=0
    ! Fill out coordinates, delta-density, and indices in a way that is
    ! reasonable fast to process.
    l=0
    do i=1,batch%n_point
        ! irreducible radial integration point
        v0=batch%radial_coordinate(:,i)
        r=norm2(v0)
        ir=1.0_r8/r
        v0=v0*ir ! scale it to norm 1


        ! delta-density at this point
        j=i+batch%semilocal_irr_offset
        dlrho=density%irr_rho(j,1)-density%irr_free_rho(j)
        ! delta-gradient projected onto radial direction? Good idea?
        dlrhoderiv=dot_product( density%irr_grad_rho(:,j,1)-density%irr_free_grad_rho(:,j),v0 )


        !dlrho=density%irr_rho(j,1)
        dlrho=density%irr_free_rho(j)
        ! delta-gradient projected onto radial direction? Good idea?
        !dlrhoderiv=dot_product( density%irr_grad_rho(:,j,1),v0 )


        ! integration weight times partition function times thing that fixes prefactor
        ispc=p%species( batch%index_atom(i) )
        irad=batch%index_radial(i)
        iang=batch%index_angular(i)
        w=grid%wts%angular_weight(ispc,irad,iang)
        w=w*batch%partition_function(i)
        w=w*weightprefactor

        ! unfold it to the full array
        do j=1,batch%unfold_ctr(i)
            jatm=batch%unfold_atom(j,i)
            ! Maybe cycle if not irreducible atom. Should speed things up.
            iop=batch%unfold_operation(j,i)
            jrad=batch%unfold_index_radial(j,i)
            ! Rotate irreducible radial integration point
            v1=matmul(sym%op(iop)%m,v0)
            ! Now we have unfolded the irreducible to the total:
            l=l+1 ! increment counter for number of points
            mih%x(l)=v1(1)
            mih%y(l)=v1(2)
            mih%z(l)=v1(3)
            mih%r(l)=r
            mih%w(l)=w
            mih%ind_atm(l)=jatm
            mih%ind_rad(l)=jrad
            mih%delta_rho(l)=dlrho
            mih%delta_rho_grad(l)=dlrhoderiv
            ! store r-value for later
            mih%rvalbuf(jrad,jatm)=r
        enddo
    enddo
    ! Make a note of the number of integration points.
    mih%n_integration_point=l
end subroutine

!> initialize the integration helper
subroutine setup_integration_helper_multipoleexpansion(mih,grid,p,sym,lmax,mw,mem)
    !> integration helper
    type(rl_integration_helper_multipoleexpansion), intent(out) :: mih
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> largest l quantum number in multipole
    integer, intent(in) :: lmax
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    integer :: l,m

    ! Set the largest angular index
    mih%n_iang=0
    do l=0,lmax
    do m=-l,l
        mih%n_iang=mih%n_iang+1
    enddo
    enddo
    ! No integration points
    mih%n_integration_point=0

    ! What is the largest i_radial?
    mih%n_radial_shell=0
    do l=1,grid%n_irr_batch
        mih%n_radial_shell=max(mih%n_radial_shell,maxval(grid%irr_batch(l)%index_radial))
    enddo
    call mw%allreduce('max',mih%n_radial_shell)

    ! Maximal number of integration points:
    l=grid%max_n_point_per_batch*sym%n_operation
    ! Space for temporary buffers
    call mem%allocate(mih%x,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%y,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%z,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%r,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%w,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%delta_rho,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%delta_rho_grad,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%Ylm,[l,mih%n_iang],scalable=.true.,persistent=.true.)
    call mem%allocate(mih%ind_atm,l,scalable=.true.,persistent=.true.)
    call mem%allocate(mih%ind_rad,l,scalable=.true.,persistent=.true.)
    mih%x=0.0_r8
    mih%y=0.0_r8
    mih%z=0.0_r8
    mih%r=0.0_r8
    mih%w=0.0_r8
    mih%delta_rho=0.0_r8
    mih%delta_rho_grad=0.0_r8
    mih%Ylm=0.0_r8
    mih%ind_atm=0
    mih%ind_rad=0

    ! And finally, space for the buffer for expansion coefficients
    call mem%allocate(mih%coeffbuf,[mih%n_iang,mih%n_radial_shell,p%n_atom],scalable=.false.,persistent=.true.)
    ! And the coordinate values
    call mem%allocate(mih%rvalbuf,[mih%n_radial_shell,p%n_atom],scalable=.false.,persistent=.true.)
    mih%coeffbuf=0.0_r8
    mih%rvalbuf=0.0_r8
end subroutine

!> destroy the helper
subroutine destroy_integration_helper_multipoleexpansion(mih,mem)
    !> integration helper
    type(rl_integration_helper_multipoleexpansion), intent(inout) :: mih
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    call mem%deallocate(mih%x,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%y,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%z,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%r,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%w,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%delta_rho,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%delta_rho_grad,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%Ylm,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%ind_atm,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%ind_rad,scalable=.true.,persistent=.true.)
    call mem%deallocate(mih%coeffbuf,scalable=.false.,persistent=.true.)
    call mem%deallocate(mih%rvalbuf,scalable=.false.,persistent=.true.)
end subroutine

!> create potential spline. Hmmm.
subroutine create_potential_spline(buf_coeff,buf_rval,grid,free,p,mw,mem)
    !> buffer with expansion coefficients
    real(r8), dimension(:,:,:), intent(in) :: buf_coeff
    !> buffer with radial coordinates
    real(r8), dimension(:,:), intent(in) :: buf_rval
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> free atom quantities
    type(rl_free_atom), intent(in) :: free
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(in) :: mem

    integer :: n_angular

    !init: block
    type(rl_hartree_multipole_spline) :: spl
    real(r8), dimension(:), allocatable :: xval,yval
    integer :: iatm,irad,iang,l,m,i,j

    ! real(r8) :: f0,f1,f2
    ! integer :: i,j,k,l
    ! integer :: iatm,irad,iang,ispc
    !
    ! f0=0.0_r8
    ! do iatm=1,1 !p%na
    !     do iang=1,1 !size(buf_coeff,1)
    !
    !     write(*,*) 'iang:',iang
    !
    !     do irad=1,size(buf_rval,1)
    !         f0=f0+buf_coeff(iang,irad,iatm)*grid%wts%radial_weight(1,irad)*buf_rval(irad,iatm)**2
    !         write(*,*) irad,buf_rval(irad,iatm),buf_coeff(iang,irad,iatm),grid%wts%radial_weight(1,irad),f0*(2*rl_twopi)
    !     enddo
    !     enddo
    ! enddo

    !init: block

        ! Number of splines per atom, sort of.
        n_angular=size(buf_coeff,1)
        allocate(xval(size(buf_coeff,2)))
        allocate(yval(size(buf_coeff,2)))


        do iatm=1,p%n_atom

            do l=0,1 !lmax
            do m=-l,l
                iang=l**2+l+1+m
                xval=0.0_r8
                yval=0.0_r8

                i=0
                do irad=1,size(buf_coeff,2)
                    if ( buf_rval(irad,iatm) .gt. rl_tiny ) then
                        i=i+1
                        xval(i)=buf_rval(irad,iatm)
                        yval(i)=buf_coeff(iang,irad,iatm)
                    endif
                enddo
                call create_spline_invloggrid_with_deriv(spl,grid,xval(1:i),yval(1:i),&
                lnum=l,i_species=p%species(iatm),free=free)

            enddo
            enddo
            do iang=1,n_angular
                ! grab coordinates and function values
                ! create a spline from this: this should be delta-density as a spline.
            enddo
        enddo
    !end block init

    ! Radially integrate the l=0 term?
    !do iatm=1,p%na
    !enddo

end subroutine

!> spline the radial delta-density thing when I know values and derivatives.
subroutine create_spline_invloggrid_with_deriv(spl,grid,x,y,lnum,i_species,free)
    !> resulting spline
    class(rl_hartree_multipole_spline), intent(out) :: spl
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> x-values
    real(r8), dimension(:), intent(in) :: x
    !> y-values
    real(r8), dimension(:), intent(in) :: y
    !> l angular index
    integer, intent(in) :: lnum
    !> species index
    integer, intent(in) :: i_species
type(rl_free_atom), intent(in) :: free

    !init: block
    integer :: i
    real(r8) :: f0,f1,f2

    !getcubicspline: block
    real(r8), dimension(:), allocatable :: xp,yp,zp,hp
    real(r8), dimension(:), allocatable :: wy,wd,wu,wl
    real(r8) :: a0,a1,a2,a3

    !real(r8) :: f0,f1,f2,r,rfac,h
    real(r8) :: r,rfac,h
    integer :: n,nx,u
    !integer :: i,j,k,l
    integer :: j,k,l

    !integratespline: block
    real(r8), dimension(:), allocatable :: belowconst,aboveconst
    real(r8), dimension(:,:), allocatable :: dr0
    real(r8) :: sra2,sra3,sra4
    !real(r8) :: a0,a1,a2,a3,b0,b1,b2,b3,b4,b5,b6,b7
    real(r8) :: b0,b1,b2,b3,b4,b5,b6,b7
    !real(r8) :: f0,f1,f2,f3,r0,r1,h,A,B
    real(r8) :: f3,r0,r1,A,B
    real(r8) :: x0,x1
    real(r8) :: integralA,integralB
    !integer :: i,j,nx

    ! Set some initial things, basic things about the spline
    !init: block
        ! Set information in the spline that we need for evaluation later
        spl%r_min  = x(1)
        spl%r_max  = maxval(x)
        spl%n_knot = size(x)
        spl%rscale = grid%wts%radial_scale(i_species)
        spl%invnsq = (1.0_r8/real(grid%wts%n_radial(i_species)+1,r8))**2
        spl%np1    = real(grid%wts%n_radial(i_species)+1,r8)
        spl%invC   = 1.0_r8/spl%rscale
        ! And space for coefficients
        allocate(spl%coeff_density(4,spl%n_knot))
        spl%coeff_density=0.0_r8
        ! Make note of which l we have
        spl%l=lnum
        ! Numerically integrate the spline to get the point-multipole coefficient.
        ! this number will also be used to normalize the spline properly later.
        spl%farfield=0.0_r8
        do i=1,spl%n_knot
            f0=y(i)*x(i)*x(i)*2*rl_twopi
            f0=f0*grid%wts%radial_weight(i_species,i)
            spl%farfield=spl%farfield+f0
        enddo
    !end block init

    ! Set up so that we can solve for the rho-spline:
    !getcubicspline: block

        nx=grid%wts%n_radial(i_species)
        n=size(x)

        ! Create padded arrays that extend slightly outside
        ! the region we originally have. Will be removed later
        ! just use it now to ensure correct limits.
        allocate(xp(0:n+1))
        allocate(yp(0:n+1))
        allocate(zp(0:n+1))
        allocate(hp(0:n+1))
        xp=0.0_r8
        yp=0.0_r8
        zp=0.0_r8
        hp=0.0_r8

        ! Convert the function into function*r^2 and spline that instead.
        do i=1,n
            f0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
            f1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )
            xp(i)=f0
            yp(i)=y(i)*f0*f0            ! function value*r^2
            hp(i)=f1-f0                 ! width of segment
        enddo

        ! fill out things for the last point:
        i=n+1
        f0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
        f1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )
        xp(i)=f0
        yp(i)=0.0_r8
        zp(i)=0.0_r8
        hp(i)=f1-f0
        ! fill out things for the zeroth point
        xp(0)=0.0_r8
        yp(0)=0.0_r8
        zp(0)=0.0_r8
        hp(0)=xp(1)

        ! Construct tridiagonal system
        allocate(wy(n+1))
        allocate(wd(n))
        allocate(wu(n-1))
        allocate(wl(n-1))
        wy=0.0_r8
        wd=0.0_r8
        wu=0.0_r8
        wl=0.0_r8

        ! First build the coefficient matrix
        do i=1,n-1
            ! Diagonal element:
            wd(i)=2*hp(i-1)+2*hp(i)
            ! superdiagonal
            wu(i)=hp(i)
            ! subdiagonal
            wl(i)=hp(i)
        enddo
        ! And the final diagonal element
        wd(n)=2*hp(n-1)+2*hp(n)
        wd=-wd/3.0_r8
        wu=-wu/3.0_r8
        wl=-wl/3.0_r8
        ! Then the vector of values to fit
        do i=1,n
            wy(i)=( yp(i)-yp(i-1) )/hp(i-1) - ( yp(i+1)-yp(i) )/hp(i)
        enddo
        ! Solve tridiagonal system
        call dgtsv(n,1,wl,wd,wu,wy(1:n),n,i)

        ! Back-substitute to get spline
        do i=1,n !-1
            a0 = yp(i)
            a1 = ( yp(i+1) - yp(i) )/hp(i) - (2*wy(i) + wy(i+1) )*hp(i)/3.0_r8
            a2 = wy(i)
            a3 = ( wy(i+1)-wy(i) )/(3*hp(i))
            spl%coeff_density(:,i)=[a0,a1,a2,a3]
        enddo
        ! For very short distances I treat it specially:
        a0=spl%coeff_density(1,1)
        a1=spl%coeff_density(2,1)
        a2=spl%coeff_density(3,1)
        a3=spl%coeff_density(4,1)
        r=spl%r_min
        spl%sa2=a2 + (6*a0 - 3*a1*r)/r**2
        spl%sa3=(-8*a0 + r*(5*a1 - 2*a2*r))/r**3
        spl%sa4=(3*a0 - 2*a1*r + a2*r**2)/r**4

        ! ! Now ... in case it is the l=0 spline, it could make sense to normalize
        ! ! it so that it analytically integrates to the correct thing. Either that
        ! ! or use the supplied radial weights to normalize it. Not sure what is the
        ! ! most stringent. Nice to check at least. First integrate the little stump:
        ! r=spl%r_min
        ! f0=(spl%sa2*r**3)/3 + (spl%sa3*r**4)/4 + (spl%sa4*r**5)/5
        ! ! Then the rest of the pieces. I have already
        ! spl%coeff_density(:,spl%n_knot)=0.0_r8
        ! do i=1,spl%n_knot
        !     a0=spl%coeff_density(1,i)
        !     a1=spl%coeff_density(2,i)
        !     a2=spl%coeff_density(3,i)
        !     a3=spl%coeff_density(4,i)
        !     f0=f0 + a0*hp(i) + (a1*hp(i)**2)/2 + (a2*hp(i)**3)/3 + (a3*hp(i)**4)/4
        ! enddo
        ! f0=f0*2*rl_twopi
        !
        ! ! Now make sure that it normalizes to what it should.
        ! f0=spl%farfield/f0
        ! spl%sa2=spl%sa2*f0
        ! spl%sa3=spl%sa3*f0
        ! spl%sa4=spl%sa4*f0
        ! spl%coeff_density=spl%coeff_density*f0
        !
        ! ! Integrate again to see if it workd.
        ! r=spl%r_min
        ! f0=(spl%sa2*r**3)/3 + (spl%sa3*r**4)/4 + (spl%sa4*r**5)/5
        ! do i=1,spl%n_knot
        !     a0=spl%coeff_density(1,i)
        !     a1=spl%coeff_density(2,i)
        !     a2=spl%coeff_density(3,i)
        !     a3=spl%coeff_density(4,i)
        !     f0=f0 + a0*hp(i) + (a1*hp(i)**2)/2 + (a2*hp(i)**3)/3 + (a3*hp(i)**4)/4
        ! enddo

        ! Dump and test.
        open(newunit=u, file='outfile_hartree_spline', status='replace', action='write')
            n=500
            r=spl%r_min*0.001_r8
            f0=spl%r_max*1.5_r8
            rfac=(f0/r)**(1.0_r8/n)
            do i=1,n
                call spl%rho_val_deriv(r,f0,f1)
                f2=free%species(1)%potential%evaluate(r)
                write(u,*) r,f0,f1,f2
                r=r*rfac
            enddo
        close(u)

    !end block getcubicspline

    ! Need to create a function so that I can evaluate the potential.
    !integratespline: block

        nx=grid%wts%n_radial(i_species)

        ! The potential is given by two integrals, one from zero to r, and one
        ! from infinity to r. It is a good idea to do this integral to the knots
        ! in the spline, since then what remains to do is just to integrate
        ! from the knots to where we currently are, much easier. There is a little
        ! algebra here to generate all the coefficients, but it's really fast.

        ! Dummy space for converted spline coefficients.
        allocate(dr0(4,spl%n_knot))
        dr0=0.0_r8

        ! First, convert the spline parameters so that they are polynomials
        ! in r, not r-r_i. Makes the rest of the algebra much much easier.
        f0=0.0_r8
        do i=1,spl%n_knot
            r0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
            r1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )
            a0=spl%coeff_density(1,i)
            a1=spl%coeff_density(2,i)
            a2=spl%coeff_density(3,i)
            a3=spl%coeff_density(4,i)

            dr0(1,i)=a0 - a1*r0 + a2*r0**2 - a3*r0**3
            dr0(2,i)=a1 - 2*a2*r0 + 3*a3*r0**2
            dr0(3,i)=a2 - 3*a3*r0
            dr0(4,i)=a3
        enddo

        ! Integrate the spline analytically and normalize, so that an analytical integration
        ! gives the same as the numerical integration to get the far-field component above.
        ! But only if we actually have a far-field component. Otherwise I have to integrate
        ! the square, or some moment or something to get the normalization right.
        if ( spl%farfield .lt. 1E-10_r8 ) then
            write(*,*) 'No far-field term. Have to consider this.'
            stop
        endif

        f0=0.0_r8
        do i=1,spl%n_knot
            r0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
            r1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )
            a0=dr0(1,i)
            a1=dr0(2,i)
            a2=dr0(3,i)
            a3=dr0(4,i)
            f0=f0-(a0*r0) - (a1*r0**2)/2 - (a2*r0**3)/3 - (a3*r0**4)/4 + a0*r1 + (a1*r1**2)/2 + (a2*r1**3)/3 + (a3*r1**4)/4
        enddo
        ! And add the little stump in the beginning to the integral.
        r0=spl%r_min
        f0=f0+(spl%sa2*r0**3)/3 + (spl%sa3*r0**4)/4 + (spl%sa4*r0**5)/5
        !f0=f0*2*rl_twopi
        f0=spl%farfield/f0
        ! Scale spline and short-range components.
        dr0=dr0*f0
        sra2=spl%sa2*f0
        sra3=spl%sa3*f0
        sra4=spl%sa4*f0

        ! Now start thinking for real. Get the constants above and below.
        allocate(belowconst(spl%n_knot))
        allocate(aboveconst(spl%n_knot))
        belowconst=0.0_r8
        aboveconst=0.0_r8
        do i=1,spl%n_knot
            a0=dr0(1,i)
            a1=dr0(2,i)
            a2=dr0(3,i)
            a3=dr0(4,i)

            x0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
            x1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )

            f0=a0*(Log(x1/x0)) + a1*(-x0 + x1) + a2*(-x0**2/2 + x1**2/2) + a3*(-x0**3/3. + x1**3/3)
            aboveconst(i)=f0

            f0=a0*(-x0 + x1)  + a1*((-x0**2 + x1**2)/2) +   a2*((-x0**3 + x1**3)/3) +   a3*((-x0**4 + x1**4)/4)
            belowconst(i)=f0
        enddo

        ! Now make space for the actual coefficients.
        allocate(spl%coeff_potential(8,spl%n_knot))
        spl%coeff_potential=0.0_r8

        ! Now start accumulating things to define a new funny spline for the potential.
        ! Not really a spline, but a piecewise function.
        do i=1,spl%n_knot
            !@TODO Add all special l's here. There are some annoying limits.

            a0=dr0(1,i)
            a1=dr0(2,i)
            a2=dr0(3,i)
            a3=dr0(4,i)
            x0=-spl%rscale*log( 1.0_r8 - (real(i,r8)/real(nx+1,r8))**2 )
            x1=-spl%rscale*log( 1.0_r8 - (real(i+1,r8)/real(nx+1,r8))**2 )

            ! Prefactor for theings that are below this segment:
            b4=0.0_r8
            do j=1,i-1
                b4=b4+belowconst(j)
            enddo
            ! And the little stump from 0 to r_min:
            b4=b4+(sra2*spl%r_min**3)/3 + (sra3*spl%r_min**4)/4 + (sra4*spl%r_min**5)/5
            ! Only works for l=0? Maybe.
            b4 = b4 + a0*x1 + (a1*x1**2)/2 + (a2*x1**3)/3 + (a3*x1**4)/4
            ! Prefactor for log(r)
            b5 = a0
            ! Normal polynomial coefficients
            b0=-a0 - a0*Log(x0) - a1*x0 - (a2*x0**2)/2 - (a3*x0**3)/3
            ! Only works for l=0
            do j=i+1,spl%n_knot-1 ! Not sure about the last one, does not count as real segment.
                b0=b0+aboveconst(j)
            enddo
            b1=a1/2
            b2=a2/6
            b3=a3/12
            ! Store coefficients
            !spl%coeff_potential(:,i)=[b0,b1,b2,b3,b4,b5,b6,b7]
            ! Evaluate for good measure.
            f0= b0 + b1*x0 + b2*x0**2 + b3*x0**3 + b4/x0 + b5*log(x0)

            write(*,*) i,x0,f0
        enddo

    !end block integratespline

write(*,*) 'continue here fit spline'
stop

end subroutine

!> evaluate radial part multipole rho expansion
elemental subroutine evaluate_rho_spline(h,r,y)
    !> function handle
    class(rl_hartree_multipole_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y

    integer :: i
    real(r8) :: a0,a1,a2,a3,alpha,r0,x

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        return
    elseif ( r .lt. h%r_min ) then
        y=h%sa2*r*r + h%sa3*r*r*r + h%sa4*r*r*r*r
        return
    endif

    ! invert the coordinates to grid index
    ! r(i) = - C ln ( 1 - (i/(n+1))**2 ), we have
    ! i(r) = (n+1) * ( 1 - exp(-r/C) )**(1/2)
    alpha = h%np1 * sqrt( 1-exp(-r*h%invC) )
    i=int(alpha)
    i=max(i,1)
    i=min(i,h%n_knot)
    r0= - h%rscale*log(1-h%invnsq*i**2)
    x=r-r0
    a0=h%coeff_density(1,i)
    a1=h%coeff_density(2,i)
    a2=h%coeff_density(3,i)
    a3=h%coeff_density(4,i)
    y = a0 + a1*x + a2*x*x + a3*x*x*x
end subroutine

!> evaluate radial part multipole rho expansion + derivative
elemental subroutine evaluate_rho_spline_deriv(h,r,y,dy)
    !> function handle
    class(rl_hartree_multipole_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y
    !> derivative
    real(r8), intent(out) :: dy

    integer :: i
    real(r8) :: a0,a1,a2,a3,alpha,r0,x

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        dy=0.0_r8
        return
    elseif ( r .le. h%r_min ) then
        y=h%sa2*r*r + h%sa3*r*r*r + h%sa4*r*r*r*r
        dy=2*h%sa2*r + 3*h%sa3*r*r + 4*h%sa4*r*r*r
        return
    endif

    ! invert the coordinates to grid index
    ! r(i) = - C ln ( 1 - (i/(n+1))**2 ), we have
    ! i(r) = (n+1) * ( 1 - exp(-r/C) )**(1/2)
    alpha = h%np1 * sqrt( 1-exp(-r*h%invC) )
    i=int(alpha)
    i=max(i,1)
    i=min(i,h%n_knot)
    r0= - h%rscale*log(1-h%invnsq*i**2)
    x=r-r0
    a0=h%coeff_density(1,i)
    a1=h%coeff_density(2,i)
    a2=h%coeff_density(3,i)
    a3=h%coeff_density(4,i)
    y = a0 + a1*x + a2*x*x + a3*x*x*x
    dy = a1 + 2*a2*x + 3*a3*x*x
end subroutine

!> evaluate radial part multipole rho expansion
elemental subroutine evaluate_potential_spline(spl,r,y)
    !> function handle
    class(rl_hartree_multipole_spline), intent(in) :: spl
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y

    integer :: i
    !real(r8) :: half=0.5_r8,third=1.0_r8/3.0_r8,fourth=0.25_r8
    real(r8) :: a0,a1,a2,a3,alpha,r0,x,h

    ! Check if we are in any if the special regions.
    if ( r .ge. spl%r_max ) then
        ! Have to think here, not sure when we can
        ! replace with the asymptotic.
        y=0.0_r8
        return
    elseif ( r .lt. spl%r_min ) then
        ! In this case it's really simple:
        y=(spl%sa2*r**2)/real(3 + spl%l,r8) + (spl%sa3*r**3)/real(4 + spl%l,r8) + (spl%sa4*r**4)/real(5 + spl%l,r8)
        return
    endif

    ! invert the coordinates to grid index
    ! r(i) = - C ln ( 1 - (i/(n+1))**2 ), we have
    ! i(r) = (n+1) * ( 1 - exp(-r/C) )**(1/2)
    ! alpha = h%np1 * sqrt( 1-exp(-r*h%invC) )
    ! i=int(alpha)
    ! i=max(i,1)
    ! i=min(i,h%n_knot)
    ! r0= - h%rscale*log(1-h%invnsq*i**2)
    ! x=r-r0
    ! a0=h%coeff_density(1,i)
    ! a1=h%coeff_density(2,i)
    ! a2=h%coeff_density(3,i)
    ! a3=h%coeff_density(4,i)
    ! y = a0 + a1*x + a2*x*x + a3*x*x*x
end subroutine



end module
