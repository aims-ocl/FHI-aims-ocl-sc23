module rlsy_free_atom
!!
!! This module handles some free-atom splined quantities
!! I just use it for debugging at this moment. Please ignore.
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_exitcode_param
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: rl_chop
use rlsy_sorting, only: rl_return_unique
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully
use rlsy_crystalstructure, only: rl_crystalstructure

implicit none
private
public :: rl_free_atom

!> radial spline for a basis function
type rl_free_atom_spline
    !> inner cutoff
    real(r8) :: r_min=-rl_huge
    !> outer cutoff
    real(r8) :: r_max=-rl_huge
    !> how many knots in the spline
    integer, private :: n_knot=-rl_hugeint
    ! pre-computed parameters to speed up evaluation:
    real(r8), private :: p1=-rl_huge   !< r_grid_inc, grid increment
    real(r8), private :: ip0=-rl_huge  !< 1/r_min
    real(r8), private :: lp0=-rl_huge  !< log(r_min)
    real(r8), private :: lp1=-rl_huge  !< log(r_grid_inc)
    real(r8), private :: ilp1=-rl_huge !< 1/log(r_grid_inc)
    !> spline coefficients (4,n_knot)
    real(r8), dimension(:,:), allocatable, private :: coeff
    contains
        ! Evaluate the free atom density at some arbitrary point.
        procedure :: evaluate => radial_function
end type

!> some free atom quantities that might come handy
type rl_free_atom_species
    !> Nuclear charge. Not sure why not integer.
    real(r8) :: Z=-rl_huge
    !> Number of electrons? Always the same as Z? Not sure. Seems likely though. Also not sure if integer.
    real(r8) :: n_electron=-rl_huge
    !> Density, defined as a spline
    type(rl_free_atom_spline) :: rho
    !> Radial charge density derivative
    type(rl_free_atom_spline) :: rho_deriv
    !> Free atom electrostatic potential.
    type(rl_free_atom_spline) :: potential
    !> Electrostatic potential at the nuclues
    real(r8) :: potential_at_nucleus=-rl_huge
    !> What is the cutoff for the free atom charge density for this species
    real(r8) :: cutoff=-rl_huge
    !> Max l in multipole expansion for this species
    integer :: l_hartree=-rl_hugeint
end type

!> Container for free atom quantities. Not all, but enough information to determine it.
type rl_free_atom
    !> how many different species are there
    integer :: n_species=-rl_hugeint
    !> free atom quantities per species
    type(rl_free_atom_species), dimension(:), allocatable :: species
    contains
        !> create object
        procedure :: generate
end type

contains

!> rearrange the information from AIMS in an easy-to-use format, just to speed up development for myself
subroutine generate(free,p,mw,mem,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero,r_grid_min,r_grid_inc,atom_radius)
    !> irreducible projection
    class(rl_free_atom), intent(out) :: free
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> mpi helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> free atom information from AIMS
    real(r8), dimension(:,:,:), intent(in) :: free_rho_spl
    real(r8), dimension(:,:,:), intent(in) :: free_drho_dr_spl
    real(r8), dimension(:,:,:), intent(in) :: free_pot_es_spl
    real(r8), dimension(:), intent(in) :: free_pot_es_at_zero
    real(r8), dimension(:), intent(in) :: r_grid_min
    real(r8), dimension(:), intent(in) :: r_grid_inc
    real(r8), dimension(:), intent(in) :: atom_radius

    !init: block
    integer, dimension(:), allocatable :: di
    integer :: ispc,iatm

    !init: block

        ! Get the number of species
        call rl_return_unique(p%species,di,mem)
        free%n_species=size(di)
        call mem%deallocate(di,persistent=.true.,scalable=.false.)

        ! Idiot check to make sure me and the rest of the world agree
        ! on how indices work. It failed once, the test remains.
        if ( maxval(p%species) .ne. free%n_species ) then
            call rl_stop_gracefully(['Clearly I do not understand indices'],rl_exitcode_param,mw%comm)
        endif

        ! Make some space
        allocate(free%species( free%n_species ))

        ! Store some really basic things
        do ispc=1,free%n_species
            ! Get nuclear charge for this species
            do iatm=1,p%n_atom
                if ( p%species(iatm) .eq. ispc ) then
                    free%species(ispc)%Z=p%atomic_number(iatm)
                    exit
                endif
            enddo
            ! Store cutoff
            free%species(ispc)%cutoff=atom_radius(ispc)
            ! Setup the splines that define the density, one analytically normalized and one verbatim copy from AIMS.
            call setup_spline(free%species(ispc)%rho,r_grid_min(ispc),r_grid_inc(ispc),atom_radius(ispc),free_rho_spl(:,:,ispc))
            call setup_spline(free%species(ispc)%rho_deriv,r_grid_min(ispc),r_grid_inc(ispc),atom_radius(ispc),free_drho_dr_spl(:,:,ispc))
            call setup_spline(free%species(ispc)%potential,r_grid_min(ispc),r_grid_inc(ispc),atom_radius(ispc),free_pot_es_spl(:,:,ispc))
            ! Make note of the free atom electrostatic potential at the nucleus. Not sure why this would be useful.
            free%species(ispc)%potential_at_nucleus=free_pot_es_at_zero(ispc)
            !@TODO l_hartree should be input
        enddo
    !end block init

    ! ! Write free atom quantities to file.
    ! dumpfreeatom: block
    !     integer, parameter :: nvals=500
    !     integer :: u
    !     integer :: i,is
    !     real(r8) :: r,rfac,f0,f1,f2,f3
    !
    !     ! Write the splines to file for inspection
    !     if ( mw%talk ) then
    !         is=1 ! First species
    !
    !         r=free%species(is)%rho%r_min*0.5_r8
    !         f0=free%species(is)%rho%r_max*1.1_r8
    !         rfac=(f0/r)**(1.0_r8/nvals)
    !         ! Make the grid a litle tighter
    !         open(newunit=u, file='outfile_free_atom_spline', status='replace', action='write')
    !             do i=1,nvals
    !                 f0=free%species(is)%rho%evaluate(r)
    !                 f1=free%species(is)%rho_deriv%evaluate(r)
    !                 f2=free%species(is)%potential%evaluate(r)
    !                 write(u,*) r,f0,f1,f2,free%species(is)%potential_at_nucleus
    !                 r=r*rfac
    !                 write(*,*) i,r
    !             enddo
    !         close(u)
    !     endif
    ! end block dumpfreeatom
end subroutine

!> evaluate free atom spline at distance r
elemental function radial_function(h,r) result(y)
    !> function handle
    class(rl_free_atom_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8) :: y

    integer :: i
    real(r8) :: a0,a1,a2,a3,alpha

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        return
    elseif ( r .lt. h%r_min ) then
        y=0.0_r8
        return
    endif

    ! Invert log grid coordinates
    !r = 1 + log(rad/h%smallest_radial_point)/log(h%radial_increment)
    alpha = 1 + log(r*h%ip0)*h%ilp1
    i = int(alpha)
    alpha=alpha-i
    a0=h%coeff(1,i)
    a1=h%coeff(2,i)
    a2=h%coeff(3,i)
    a3=h%coeff(4,i)
    y = a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha

    ! Derivative not too tricky.
end function

! !> analytically integrate radial spline*r^2 over all of space. Well not all-all of space, but you get the point. Also 4*pi missing for some reason.
! pure function integrate_spline_times_rsq(h) result(y)
!     !> spline handle
!     class(rl_free_atom_spline), intent(in) :: h
!     !> integral
!     real(r8) :: y
!
!     real(r8) :: p0,p1,p1c,pr0,a0,a1,a2,a3,f1
!     real(r8) :: tr0,tr1,tr2
!     real(r8) :: lp1,lp2,lp3
!     integer :: i
!
!     ! pre-compute some common factors in the analytical integration
!     p0=h%smallest_radial_point
!     p1=h%radial_increment
!     p1c=p1**3
!     lp1=log(p1)
!     lp2=lp1**2
!     lp3=lp1**3
!     pr0=(p0**3)/(27*lp3)
!     y=0.0_r8
!     do i=1,h%n_pt
!         ! fetch the spline coefficients
!         a0=h%spline_coefficients(1,i)
!         a1=h%spline_coefficients(2,i)
!         a2=h%spline_coefficients(3,i)
!         a3=h%spline_coefficients(4,i)
!         ! Integrate this spline segment analytically.
!         tr0=pr0*p1**(3*i-3)
!         tr1=2*a3 - 2*a2*lp1 + 3*a1*lp2 - 9*a0*lp3
!         tr2=-2*a3 + 2*(a2 + 3*a3)*lp1 - 3*(a1 + 2*a2 + 3*a3)*lp2 + 9*(a0 + a1 + a2 + a3)*lp3
!         f1=tr0*( tr1 + tr2*p1c )
!         y=y+f1
!     enddo
! end function

! !> evaluate the radial basis functions given a 3D point
! elemental function radial_function_spline_scalar_invr(h,rad) result(y)
!     !> function handle
!     class(rl_lcao_basis_set_spline), intent(in) :: h
!     !> distance of point
!     real(r8), intent(in) :: rad
!     !> function value
!     real(r8) :: y
!
!     integer :: i
!     real(r8) :: r,a0,a1,a2,a3
!
!     ! Some safety checks here, not sure what they mean though
!     ! Function should be nothing when distance too long
!     ! Also probably zero when too short.
!     if ( rad .ge. h%cutoff ) then
!         y=0.0_r8
!         return
!     elseif ( rad .lt. h%smallest_radial_point*1E-10_r8 ) then
!         y=0.0_r8
!         return
!     endif
!
!     ! Invert log grid coordinates
!     r = 1 + log(rad/h%smallest_radial_point)/log(h%radial_increment)
!     i = int(r)
!     r=r-i
!     a0=h%spline_coefficients(1,i)
!     a1=h%spline_coefficients(2,i)
!     a2=h%spline_coefficients(3,i)
!     a3=h%spline_coefficients(4,i)
!     y = a0 + a1*r + a2*r*r + a3*r*r*r
!     y = y/rad
! end function

!> get the spline from AIMS to a more convenient form, defined everywhere and all that, with smooth derivatives.
subroutine setup_spline(spline,r_grid_min,r_grid_inc,cutoff,coeff)
    !> resulting spline
    type(rl_free_atom_spline), intent(out) :: spline
    !> smallest grid point
    real(r8), intent(in) :: r_grid_min
    !> increment to build grid
    real(r8), intent(in) :: r_grid_inc
    !> cutoff for the spline
    real(r8), intent(in) :: cutoff
    !> spline coefficients
    real(r8), dimension(:,:), intent(in) :: coeff

    !donottaper: block
        real(r8) :: r
        integer :: i

    !init: block
        ! Store some basic things about the spline, and pre-computed
        ! constants for faster evaluation later.
        spline%r_max=cutoff
        spline%r_min=r_grid_min
        spline%p1=r_grid_inc
        spline%ip0=1.0_r8/spline%r_min
        spline%lp0=log(spline%r_min)
        spline%lp1=log(spline%p1)
        spline%ilp1=1.0_r8/spline%lp1
    !end block init

    !donottaper: block
        ! Just copy the spline instead.

        ! Copy only the necessary spline coefficients.
        spline%n_knot=0
        do i=1,size(coeff,2)
            r=exp( (i-1)*spline%lp1+spline%lp0 )
            if ( r .gt. cutoff+1E-5_r8 ) then
                exit
            else
                spline%n_knot=spline%n_knot+1
            endif
        enddo
        allocate(spline%coeff(4,spline%n_knot))
        spline%coeff=coeff(:,1:spline%n_knot)
    !end block donottaper
end subroutine

end module
