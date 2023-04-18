module rlsy_basis_set
!!
!! This module handles basis functions: both in terms of evaluating and organizing them.
!! Also, the routines that generate rotation matrices are stored here.
!! This will likely be stripped down to the bare minimum at some point.
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_pi,rl_twopi,rl_sqtol,rl_exitcode_param,rl_exitcode_symmetry,&
                          rl_exitcode_blaslapack,rl_exitcode_memory
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: rl_chop,rl_determ,tochar
use rlsy_sorting, only: rl_return_unique,rl_qsort
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup_operation

implicit none
private
public :: rl_lcao_basis_set

!> information about the basis set attached to an atom. Some information is redundant for easy of access and shorter code elsewhere.
type rl_lcao_basis_set_species
    !> how many basis functions live here
    integer :: n_basis=-rl_hugeint
    !> radial index
    integer, dimension(:), allocatable :: basis_fn
    !> l-number
    integer, dimension(:), allocatable :: basis_l
    !> m-number
    integer, dimension(:), allocatable :: basis_m
    !> combined l and m number ( = 1+l+l^2+m ) in case you only want to lookup once!
    integer, dimension(:), allocatable :: basis_ang
    !> how many groups of constant l are there
    integer :: n_basis_group=-rl_hugeint
    !> l quantum number for each group (n_basis_groups)
    integer, dimension(:), allocatable :: basis_group_l
    !> groups of basis functions (2,n_basis_groups_on_atom), start and end of indices
    integer, dimension(:,:), allocatable :: basis_group
    !> radial cutoff per basis function. Same information as in radial below, but accessible without a million indices.
    real(r8), dimension(:), allocatable :: basis_cutoff
end type

!> information about the basis set attached to an atom. Some information is redundant for easy of access and shorter code elsewhere.
type rl_lcao_basis_set_atom
    !> how many basis functions live here
    integer :: n_basis_on_atom=-rl_hugeint
    !> what is the offset from this atom? with respect to the unit cell, that is.
    integer :: basis_index_offset=-rl_hugeint
    !> which basis functions live here, in terms of the global index in the unit cell. Not sure if used.
    integer, dimension(:), allocatable :: basis_on_atom
    !> radial index
    integer, dimension(:), allocatable :: basis_fn
    !> l-number
    integer, dimension(:), allocatable :: basis_l
    !> m-number
    integer, dimension(:), allocatable :: basis_m
    !> combined l and m number ( = 1+l+l^2+m ) in case you only want to lookup once!
    integer, dimension(:), allocatable :: basis_ang
    !> how many groups of constant l are there
    integer :: n_basis_groups_on_atom=-rl_hugeint
    !> l quantum number for each group (n_basis_groups_on_atom)
    integer, dimension(:), allocatable :: basis_group_l
    !> groups of basis functions (2,n_basis_groups_on_atom), start and end of indices
    integer, dimension(:,:), allocatable :: basis_group
    !> radial cutoff per basis function. Same information as in radial below, but accessible without a million indices.
    real(r8), dimension(:), allocatable :: basis_cutoff
end type

!> radial spline for a basis function
type rl_lcao_basis_set_spline
    !> which species does this basis function belong to?
    integer :: i_species=-rl_hugeint
    !> cutoff in large cartesian distance, Cartesian coordinates
    real(r8) :: r_max=-rl_huge
    !> cutoff where we switch from spline to polynomial, Cartesian coordinates
    real(r8) :: r_min=-rl_huge
    !> how many knots in the spline
    integer, private :: n_knot=-rl_hugeint
    ! pre-computed parameters to speed up evaluation:
    real(r8), private :: p1=-rl_huge   !< r_grid_inc, grid increment
    real(r8), private :: ip0=-rl_huge  !< 1/r_min
    real(r8), private :: lp0=-rl_huge  !< log(r_min)
    real(r8), private :: lp1=-rl_huge  !< log(r_grid_inc)
    real(r8), private :: ilp1=-rl_huge !< 1/log(r_grid_inc)
    real(r8), private :: lb1=-rl_huge  !< coefficient for small r
    real(r8), private :: lb2=-rl_huge  !< coefficient for small r^2
    real(r8), private :: lb3=-rl_huge  !< coefficient for small r^3
    real(r8), private :: lb4=-rl_huge  !< coefficient for small r^4
    !> spline coefficients (4,n_knot)
    real(r8), dimension(:,:), allocatable, private :: coeff
    real(r8), dimension(:,:), allocatable, private :: coeff_deriv
    contains
        !> get value at arbitrary r
        procedure :: val => radial_function
        !> get value+derivative at arbitrary r
        procedure :: val_der => radial_function_derivative
        !> get value+derivative+hessian at arbitrary r
        procedure :: val_der_hess => radial_function_derivative_hessian
end type

!> contains information how to project basis functions
type rl_lcao_basis_set
    !> number of basis functions per unit cell
    integer :: n_basis=-rl_hugeint
    !> number of different radial basis functions
    integer :: n_radial=-rl_hugeint
    !> number of different angular basis functions
    integer :: n_angular=-rl_hugeint
    !> highest l angular number
    integer :: l_max=-rl_hugeint
    !> basis set information per species
    type(rl_lcao_basis_set_species), dimension(:), allocatable :: species
    !> radial basis functions, guaranteed to be smooth and all that.
    type(rl_lcao_basis_set_spline), dimension(:), allocatable :: radial
    !> radial basis functions, verbatim copy from AIMS
    type(rl_lcao_basis_set_spline), dimension(:), allocatable :: radial_raw
    !> longest confinement cutoff radius in the basis set
    real(r8) :: longest_cutoff=-rl_huge
    !> offset in index per atom in the unit cell (p%n_atom)
    integer, dimension(:), allocatable :: offset
    contains
        !> create object
        procedure :: generate
        !> rotation matrix, realspace
        procedure :: realspace_rotation_matrix_for_one_atom
        !> rotation matrix, reciprocal space
        procedure :: kspace_rotation_matrix
        !> rotation matrix for a block of constant l
        procedure, nopass :: rotation_matrix_constant_l => olles_rotmat
        !> evaluate spherical harmonics, derivatives and hessians.
        procedure, nopass :: Ylm => spherical_harmonic
        procedure, nopass :: Ylm_grad => spherical_harmonic_gradient
        procedure, nopass :: Ylm_grad_hess => spherical_harmonic_gradient_hessian
        !> measure size in memory, approximately
        procedure :: size_in_mem
        !> destroy object
        procedure :: destroy
end type

contains

!> take information about the basis set from AIMS and store it in a convenient form.
! eventually I might retire this, but makes development so much faster for now.
subroutine generate(basis,p,mw,mem,&
    basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
    r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl)
    !> irreducible projection
    class(rl_lcao_basis_set), intent(out) :: basis
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> mpi helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> basis set information from AIMS
    integer, dimension(:), intent(in) :: basis_m
    integer, dimension(:), intent(in) :: basis_l
    integer, dimension(:), intent(in) :: basis_fn
    integer, dimension(:), intent(in) :: basis_atom
    real(r8), dimension(:), intent(in) :: outer_radius
    real(r8), dimension(:), intent(in) :: r_grid_min
    real(r8), dimension(:), intent(in) :: r_grid_inc
    real(r8), dimension(:,:,:), intent(in) :: basis_wave_spl
    real(r8), dimension(:,:,:), intent(in) :: basis_deriv_spl

    !init: block
    integer, dimension(:), allocatable :: stoa,di
    integer :: s1
    integer :: a1,i,j,k

    !sortintogroups: block
    integer, dimension(:,:), allocatable :: dj1,dj2
    integer, dimension(:), allocatable :: di1
    !integer :: s1,i,j,k,l
    integer :: l

    !sortradial: block
    !integer :: i,j,k

    !sortangular: block
    !    integer :: i,l,m
    integer :: m

    ! I have to sort the basis functions a little, put them in the right place and all that.
    ! I have to sort them by atom, then group them by constant l, and store the start
    ! and stop index for each group in terms of the global index. This is information needed to
    ! construct the transformation matrices for the basis functions. This should be take almost no
    ! time and memory.
    !init: block

        ! number of basis function per unit cell
        basis%n_basis=size(basis_l)

        ! Space per species to store information
        allocate(basis%species(p%n_species))
        ! Space for offset
        allocate(basis%offset(p%n_atom))
        basis%offset=0
        ! Get an atom from species index
        call mem%allocate(stoa,p%n_species,persistent=.false.,scalable=.false.)
        do s1=1,p%n_species
            do a1=1,p%n_atom
                if ( p%species(a1) .eq. s1 ) then
                    stoa(s1)=a1
                    exit
                endif
            enddo
        enddo

        call mem%allocate(di,p%n_atom,persistent=.false.,scalable=.false.)
        di=0
        ! count basis functions per atom
        do i=1,basis%n_basis
            a1=basis_atom(i)
            di(a1)=di(a1)+1
        enddo

        ! make space for basis specifications, reset counter
        do s1=1,p%n_species
            basis%species(s1)%n_basis=di( stoa(s1) )
            allocate(basis%species(s1)%basis_fn    ( basis%species(s1)%n_basis ))
            allocate(basis%species(s1)%basis_l     ( basis%species(s1)%n_basis ))
            allocate(basis%species(s1)%basis_m     ( basis%species(s1)%n_basis ))
            allocate(basis%species(s1)%basis_cutoff( basis%species(s1)%n_basis ))
            allocate(basis%species(s1)%basis_ang   ( basis%species(s1)%n_basis ))
            basis%species(s1)%basis_fn=0
            basis%species(s1)%basis_l=0
            basis%species(s1)%basis_m=0
            basis%species(s1)%basis_cutoff=0.0_r8
            basis%species(s1)%basis_ang=0
            basis%species(s1)%n_basis=0
        enddo
        ! store specification
        di=0
        do i=1,basis%n_basis
            a1=basis_atom(i)
            di(a1)=di(a1)+1
            do s1=1,p%n_species
                if ( a1 .eq. stoa(s1) ) then
                    basis%species(s1)%n_basis=basis%species(s1)%n_basis+1
                    j=basis%species(s1)%n_basis
                    basis%species(s1)%basis_fn(j) = basis_fn(i)
                    basis%species(s1)%basis_l(j)  = basis_l(i)
                    basis%species(s1)%basis_m(j)  = basis_m (i)
                    basis%species(s1)%basis_ang(j)  = 1 + basis_l(i) + basis_l(i)**2 + basis_m(i)
                    basis%species(s1)%basis_cutoff(j)=outer_radius( basis_fn(i) )
                endif
            enddo
        enddo

        ! Get the offset per atom
        j=0
        do a1=1,p%n_atom
            s1=p%species(a1)
            basis%offset(a1)=j
            j=j+basis%species(s1)%n_basis
        enddo

        ! Now double-check that this makes sense:
        k=0
        do a1=1,p%n_atom
            s1=p%species(a1)
            do i=1,basis%species(s1)%n_basis
                j=basis%offset(a1)+i ! global index
                if ( basis_fn(j) .ne. basis%species(s1)%basis_fn(i) ) k=k+1
                if ( basis_l(j)  .ne. basis%species(s1)%basis_l(i)  ) k=k+1
                if ( basis_m(j)  .ne. basis%species(s1)%basis_m(i)  ) k=k+1
            enddo
        enddo
        if ( k .ne. 0 ) then
            call rl_stop_gracefully(['Could not index basis properly'],rl_exitcode_symmetry,mw%comm)
        endif

        call mem%deallocate(stoa,persistent=.false.,scalable=.false.)
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    !end block init

    ! On each species, I need to group the basis functions by constant
    ! radial and l, since those are the things that can transform
    ! to each other. It is quite easy, since AIMS returns the functions
    ! ordered this way to start with, always. I hope. If not I need to
    ! revise this thoroughly.
    !sortintogroups: block

        do s1=1,p%n_species
            call mem%allocate(di1,basis%species(s1)%n_basis,persistent=.false.,scalable=.false.)
            call mem%allocate(dj1,[2,basis%species(s1)%n_basis],persistent=.false.,scalable=.false.)
            di1=0
            dj1=0
            ! Get the unique combinations of radial and l index
            do i=1,basis%species(s1)%n_basis
                dj1(:,i)=[  basis%species(s1)%basis_fn(i), basis%species(s1)%basis_l(i)]
            enddo
            call rl_return_unique(dj1,dj2,mem)

            basis%species(s1)%n_basis_group=size(dj2,2)
            allocate( basis%species(s1)%basis_group(2,basis%species(s1)%n_basis_group) )
            allocate( basis%species(s1)%basis_group_l(basis%species(s1)%n_basis_group) )
            basis%species(s1)%basis_group=0
            basis%species(s1)%basis_group_l=0

            ! Get start and stop index of each group
            di1=0
            il1: do i=1,basis%species(s1)%n_basis
                do j=1,size(dj2,2)
                    if ( sum(abs(dj1(:,i)-dj2(:,j))) .eq. 0 ) then
                        di1(i)=j
                        cycle il1
                    endif
                enddo
            enddo il1
            basis%species(s1)%basis_group(1,1)=1
            do i=1,basis%species(s1)%n_basis-1
                if ( di1(i) .ne. di1(i+1) ) then
                    ! new group
                    basis%species(s1)%basis_group(2,di1(i))=i
                    basis%species(s1)%basis_group(1,di1(i)+1)=i+1
                endif
            enddo
            basis%species(s1)%basis_group(2,basis%species(s1)%n_basis_group)=basis%species(s1)%n_basis

            ! Stow away l quantum number per group
            do i=1,basis%species(s1)%n_basis_group
                basis%species(s1)%basis_group_l(i)=dj2(2,i)
            enddo

            ! Check that I got it right
            do i=1,basis%species(s1)%n_basis_group
                k=0
                do j=basis%species(s1)%basis_group(1,i),basis%species(s1)%basis_group(2,i)
                    k=k+1
                enddo
                l=basis%species(s1)%basis_group_l(i)
                if ( k .ne. 2*l+1 ) then
                    call rl_stop_gracefully(['Clearly I can not count'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo

            call mem%deallocate(di1,persistent=.false.,scalable=.false.)
            call mem%deallocate(dj1,persistent=.false.,scalable=.false.)
            call mem%deallocate(dj2,persistent=.true.,scalable=.false.)
        enddo

    !end block sortintogroups

    ! While I'm at it, could be a good idea to copy information about the radial part.
    ! Will likely be removed later once I understand things better, but for now it's nice
    ! to have it. Also a much easier way to evaluate the basis functions as opposed to
    ! going 50 subroutines deep to do so in AIMS.
    !sortradial: block

        ! this is the number of different radial functions
        basis%n_radial = maxval(basis_fn)
        allocate(basis%radial( basis%n_radial ))
        allocate(basis%radial_raw( basis%n_radial ))

        radl: do i=1,basis%n_radial
            ! Figure out which species this guy is associated with.
            do j=1,basis%n_basis
                if ( basis_fn(j) .eq. i ) then
                    k=p%species(basis_atom(j))
                    ! This is a clean spline I will orthonormalize analytically later, for debugging and such
                    call setup_spline(basis%radial(i),r_grid_min(k),r_grid_inc(k),outer_radius(i),basis_wave_spl(:,:,i),basis_deriv_spl(:,:,i),k,taper=.false.)
                    ! Verbatim copy of the splines from AIMS.
                    call setup_spline(basis%radial_raw(i),r_grid_min(k),r_grid_inc(k),outer_radius(i),basis_wave_spl(:,:,i),basis_deriv_spl(:,:,i),k,taper=.false.)
                    ! And we are done with this radial function
                    cycle radl
                endif
            enddo
        enddo radl
        ! Also store the longest radial cutoff in the system
        basis%longest_cutoff = maxval(outer_radius)
    !end block sortradial

    ! It is nice with a flat indexing of angular functions, useful sometimes
    !sortangular: block
        !integer :: i,l,m

        ! Do it per species:
        basis%l_max=0
        do i=1,p%n_species
            basis%l_max=max(basis%l_max,maxval(basis%species(i)%basis_l))
        enddo
        i=0
        do l=0,basis%l_max
        do m=-l,l
            i=i+1
        enddo
        enddo
        basis%n_angular=i
    !end block sortangular

    !     ! Test derivatives. Seems to work fine. Hmm.
    !     testderiv: block
    !         real(r8), dimension(100) :: dumval
    !         real(r8) :: sc(11,3)
    !         real(r8) :: sv(11),sr(11)
    !         real(r8), dimension(3) :: v0,v1
    !         real(r8), dimension(3) :: dr1,dr2
    !         real(r8) :: xh,yh,zh,r,ir
    !         real(r8) :: f0,f1,f2
    !         integer :: indr,inda
    !         integer :: i,ix,il,im,iang,ctr
    !
    !         !v0=[0.1_r8,0.2_r8,0.3_r8]+2.0_r8
    !         !v0=[1,2,3]
    !         call random_number(v0)
    !         xh=v0(1)/norm2(v0)
    !         yh=v0(2)/norm2(v0)
    !         zh=v0(3)/norm2(v0)
    !         r=norm2(v0)
    !         ir=1.0_r8/norm2(v0)
    !
    ! !         ! Test radial first
    ! !         indr=6
    ! !         sv=0.0_r8
    ! !         sc=0.0_r8
    ! !         do ix=1,3
    ! !             call lo_centraldifference(5,v0(ix),1E-6_r8,sc)
    ! !             do i=1,11
    ! !                 v1=v0
    ! !                 v1(ix)=sc(i,1)
    ! !                 call basis%radial(indr)%val(norm2(v1),sv(i))
    ! !             enddo
    ! !             dr1(ix)=sum(sv*sc(:,2))
    ! !         enddo
    ! !         ! Analytical derivative
    ! !         call basis%radial(indr)%val_der(norm2(v1),f0,f1)
    ! !
    ! !
    ! ! write(*,*) 'dr1',dr1,norm2(dr1)
    ! ! write(*,*) 'dr2',v0*f1/norm2(v0),f1
    ! !
    ! !
    ! !
    ! !         write(*,*) 'dr1',dr1,norm2(dr1)
    ! !         write(*,*) 'dr2',dr2,norm2(dr2)
    ! !
    ! !         ! Combined guy. Hmm.
    ! !         do ix=1,3
    ! !             call lo_centraldifference(5,v0(ix),1E-6_r8,sc)
    ! !             do i=1,11
    ! !                 v1=v0
    ! !                 v1(ix)=sc(i,1)
    ! !                 xh=v1(1)/norm2(v1)
    ! !                 yh=v1(2)/norm2(v1)
    ! !                 zh=v1(3)/norm2(v1)
    ! !                 r=norm2(v1)
    ! !                 ir=1.0_r8/norm2(v1)
    ! !                 call basis%Ylm(xh,yh,zh,inda,sv(i))
    ! !                 call basis%radial(indr)%val(r,sr(i))
    ! !             enddo
    ! !             dr1(ix)=sum(sv*sr*sc(:,2))
    ! !         enddo
    ! !         v1=v0
    ! !         xh=v1(1)/norm2(v1)
    ! !         yh=v1(2)/norm2(v1)
    ! !         zh=v1(3)/norm2(v1)
    ! !         r=norm2(v1)
    ! !         ir=1.0_r8/norm2(v1)
    ! !         call basis%Ylm_grad(xh,yh,zh,ir,inda,f0,dr2(1),dr2(2),dr2(3))
    ! !         call basis%radial(indr)%val_der(r,f1,f2)
    ! !
    ! !         write(*,*) 'dr1:',dr1
    ! !         v0=f1*dr2+f0*v1*ir*f2
    ! !         write(*,*) 'dr2:',v0
    !
    !         ! Test the angular, agressively.
    !
    ! do inda=1,9
    !         ! Now go for the angular.
    !         !inda=4
    !         do ix=1,3
    !             call lo_centraldifference(5,v0(ix),1E-4_r8,sc)
    !             do i=1,11
    !                 v1=v0
    !                 v1(ix)=sc(i,1)
    !                 xh=v1(1)/norm2(v1)
    !                 yh=v1(2)/norm2(v1)
    !                 zh=v1(3)/norm2(v1)
    !                 r=norm2(v1)
    !                 ir=1.0_r8/norm2(v1)
    !                 call basis%Ylm(xh,yh,zh,inda,sv(i))
    !                 !call basis%radial(indr)%val(norm2(v1),sv(i))
    !             enddo
    !             dr1(ix)=sum(sv*sc(:,2))
    !         enddo
    !         ! Analytical derivative
    !         v1=v0
    !         xh=v1(1)/norm2(v1)
    !         yh=v1(2)/norm2(v1)
    !         zh=v1(3)/norm2(v1)
    !         r=norm2(v1)
    !         ir=1.0_r8/norm2(v1)
    !         call basis%Ylm_grad(xh,yh,zh,ir,inda,f0,dr2(1),dr2(2),dr2(3))
    !
    ! write(*,*) ''
    ! write(*,*) 'inda:',inda
    !         write(*,*) 'dr1',dr1,norm2(dr1)
    !         write(*,*) 'dr2',dr2,norm2(dr2)
    !         write(*,*) dr1-dr2
    !
    ! enddo
    !
    ! call random_number(v0)
    ! xh=v0(1)/norm2(v0)
    ! yh=v0(2)/norm2(v0)
    ! zh=v0(3)/norm2(v0)
    ! r=norm2(v0)
    ! ir=1.0_r8/norm2(v0)
    !
    ! write(*,*) ''
    ! write(*,*) 'TEST WRT TO AIMS'
    ! write(*,*) ''
    ! ctr=0
    ! dumval=0.0_r8
    ! do il=0,6
    ! do im=-il,il
    !     ctr=ctr+1
    !
    ! write(*,*) ''
    ! write(*,*) 'point',il,im
    !     call aims_spherical_harm(v0,il,im,f0,dr2(1),dr2(2),dr2(3))
    !     iang=il**2+il+1+im
    !     call basis%Ylm_grad(xh,yh,zh,ir,iang,f1,dr1(1),dr1(2),dr1(3))
    !
    ! write(*,*) f0,dr1,norm2(dr1)
    ! write(*,*) f1,dr2,norm2(dr2)
    !
    !     dumval(ctr)=sum(abs(dr1-dr2))
    !
    ! enddo
    ! enddo
    !
    ! do i=1,ctr
    !     write(*,*) i,dumval(i)
    ! enddo
    !
    ! write(*,*) 'DONE HERE TEST DERIVATIVE'
    ! call mw%destroy()
    ! stop
    !
    !     end block testderiv

    ! Check that memory is cleared properly.
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
end subroutine

!> get the spline from AIMS to a more convenient form, defined everywhere and all that, with smooth derivatives.
subroutine setup_spline(spline,r_grid_min,r_grid_inc,cutoff,coeff,coeff_deriv,i_species,taper)
    !> resulting spline
    type(rl_lcao_basis_set_spline), intent(out) :: spline
    !> smallest grid point
    real(r8), intent(in) :: r_grid_min
    !> increment to build grid
    real(r8), intent(in) :: r_grid_inc
    !> cutoff for the spline
    real(r8), intent(in) :: cutoff
    !> spline coefficients
    real(r8), dimension(:,:), intent(in) :: coeff
    !> spline coefficients for derivatives
    real(r8), dimension(:,:), intent(in) :: coeff_deriv
    !> which species is the basis function attached to
    integer, intent(in) :: i_species
    !> make sure the function tapers off towards the cutoff?
    logical, intent(in) :: taper

    !smallr: block
    integer, parameter :: nsp=4         ! Order of spline
    integer, parameter :: ncon=1        ! Number of constraints
    integer, parameter :: niter=4       ! Number of iterations in the fitting procedure
    real(r8), dimension(ncon,nsp) :: wC ! Space for constraints
    real(r8), dimension(ncon) :: wD1     ! Space for constraints
    real(r8), dimension(nsp) :: wsol    ! Space for solution
    real(r8), dimension(:,:), allocatable :: wA     ! Coefficient for polynomial
    real(r8), dimension(:), allocatable :: wB       ! y-values to fit
    real(r8), dimension(:), allocatable :: work     ! workspace for solver
    real(r8) :: r,s_work(1),f0,f1
    integer :: n,i,iter,i_pmax,n_pbuf,lwork

    !taperspline: block
    integer, parameter :: ntaper=1      ! How many kinds of tapering will I try?
    integer, parameter :: lowrbuf=4     ! Buffer to smooth to the polynomial at low r
    real(r8), dimension(:,:), allocatable :: cleancoeff,yt !,wA,wC,wB
    real(r8), dimension(:), allocatable :: x,y,wy,wd,wu,wl
    !real(r8) :: r,f0,f1,a0,a1,a2,a3,b,c,rc,y0,yp0,ypp0
    real(r8) :: a0,a1,a2,a3,b,c,rc,y0,yp0,ypp0
    !integer :: i,l,n,m,irc,iter
    integer :: l,m,irc

    !donottaper: block
    !real(r8) :: r
    !integer :: i

    !init: block
        ! Store some basic things about the spline, and pre-computed
        ! constants for faster evaluation later.
        spline%i_species=i_species
        spline%r_max=cutoff
        spline%r_min=r_grid_min
        spline%p1=r_grid_inc
        spline%ip0=1.0_r8/spline%r_min
        spline%lp0=log(spline%r_min)
        spline%lp1=log(spline%p1)
        spline%ilp1=1.0_r8/spline%lp1
    !end block init

    ! First we try to take care of the behavior at very small r, close to
    ! the cutoff.
    !smallr: block

        ! Starting parameters for the fit.
        i_pmax=0  ! Smallest amount of points needed.
        n_pbuf=4  ! buffer above max to include in fit.

        do iter=1,niter
            ! First, get the current number of points included
            i_pmax=i_pmax+1
            n_pbuf=n_pbuf+1
            n=i_pmax+n_pbuf

            ! Reset workspace
            allocate(wA(n,nsp))
            allocate(wB(n))
            wA=0.0_r8
            wB=0.0_r8
            wC=0.0_r8
            wD1=0.0_r8
            wsol=0.0_r8

            ! Now build the design matrix
            do i=1,n
                ! x-value at this point
                r=spline%r_min*spline%p1**(i-1)
                wA(i,:)=[r,r**2,r**3,r**4]
                ! y-value at this point
                wB(i)=coeff(1,i)
            enddo

            ! And construct the constraints
            r=spline%r_min*spline%p1**(i_pmax-1)
            wC(1,:)=[r,r**2,r**3,r**4]
            wD1(1)=coeff(1,i_pmax)

            ! Solve constrained lsq problem for polynomial coefficients:
            lwork=-1
            call dgglse(n,nsp,1,wA,n,wC,1,wB,wD1,wsol,s_work,lwork,i)
            if ( i .eq. 0 ) then
                lwork = int(anint(s_work(1)))
                allocate(work(lwork))
                call dgglse(n,nsp,1,wA,n,wC,1,wB,wD1,wsol,work,lwork,i)
                deallocate(work)
            endif
            if ( i .ne. 0 ) then
                call rl_stop_gracefully(['dgglse exit status '//tochar(i)],rl_exitcode_blaslapack)
            endif

            ! Now measure the error? Think I will go with rms error
            f0=0.0_r8
            do i=1,n
                r=spline%r_min*spline%p1**(i-1)
                f1=wsol(1)*r + wsol(2)*r**2 + wsol(3)*r**3 + wsol(4)*r**4
                f0=f0+( (f1-coeff(1,i))/coeff(1,i) )**2
            enddo
            f0=sqrt(f0/real(n,r8))

            ! This is a sensibly tight tolerance, 1E-8 in relative error.
            if ( f0 .gt. 1E-8_r8 .or. iter .eq. niter ) then
                exit
            endif

            deallocate(wA)
            deallocate(wB)
        enddo
        ! That might have been somewhat reasonable? Yes no maybe.
        ! Store the polynomial coefficients:
        spline%lb1=wsol(1)
        spline%lb2=wsol(2)
        spline%lb3=wsol(3)
        spline%lb4=wsol(4)
    !end block smallr

    ! First, adjust the input spline such that it goes to zero
    ! properly at the cutoff. I do this by tapering it with a
    ! cos(x) from some arbitrary onset. Easiest way to do it is
    ! just to taper the function values and re-fit the spline.
    if ( taper ) then
    !taperspline: block

        n=size(coeff,2) ! Input number of points
        m=n+lowrbuf     ! Buffer at low r, just for creating the spline.

        ! check that the spline is sensible
        if ( n .le. 10 ) then
            call rl_stop_gracefully(['10 or less spline coefficients on input, can not be good.'],rl_exitcode_param)
        endif

        ! Some temporary workspace when tapering the function off.
        allocate(x(m))
        allocate(y(m))
        allocate(yt(n,ntaper))
        allocate(cleancoeff(4,m))
        allocate(wy(m))
        allocate(wd(m))
        allocate(wu(m-1))
        allocate(wl(m-1))
        cleancoeff=0.0_r8
        x=0.0_r8
        y=0.0_r8
        yt=0.0_r8
        wy=0.0_r8
        wd=0.0_r8
        wu=0.0_r8
        wl=0.0_r8

        ! Get parameters at the cutoff radius to help with tapering.
        r = 1 + log(spline%r_max*spline%ip0)*spline%ilp1
        i = int(r)
        r=r-i
        if ( abs(r) .gt. 1E-10_r8 ) then
            call rl_stop_gracefully(['Cutoff not at knot, I expected that.'],rl_exitcode_param)
        else
            irc=i ! Index in spline for the cutoff.
        endif
        a0=coeff(1,irc)
        a1=coeff(2,irc)
        a2=coeff(3,irc)
        a3=coeff(4,irc)
        rc=spline%r_max
        y0=a0                                                      ! y at rc
        yp0=a1*spline%ilp1/rc                                      ! y' at rc
        ypp0=(2*a2-a1*spline%lp1 )*spline%ilp1**2/spline%r_max**2  ! y'' at rc
        !yppp0=(6*a3-6*a2*spline%lp1+2*a1*spline%lp1**2)*spline%ilp1**3/r0**3 ! y''' at rsoft

        ! ideally the value and the first two derivatives should be zero at rc.
        ! Try several ways of tapering off the function. eventually. Just one for now.
        do iter=1,ntaper
            select case(iter)
            case(1)
                ! Taper with a*r*exp(-(r-b)^2/c^2), seems sensible.
                f0=y0**2 + rc**2*yp0**2 - rc**2*y0*ypp0
                if ( f0 .lt. 1E-20_r8 ) then
                    call rl_stop_gracefully(['Please revisit tapering function.'],rl_exitcode_param)
                endif
                ! Parameters for the tapering function
                b=(rc*(2*f0 + rc*(yp0*(y0 - rc*yp0) + rc*y0*ypp0)))/f0
                c=(sqrt(2.0_r8)*rc*y0)/sqrt(f0)
                ! Evaluate tapering function
                f0=exp( (rc - b)**2/c**2)
                do i=1,irc
                    r=spline%r_min*spline%p1**(i-1)
                    !r=exp( (i-1)*spline%lp1+spline%lp0 )
                    if ( abs(r-b) .lt. 10*c ) then
                        yt(i,iter)=(f0*r*y0/rc)*exp(-(r-b)**2/c**2)
                    endif
                enddo
            end select
        enddo

        ! Add y-values from polynomial at really small r.
        l=0
        do i=-lowrbuf+1,lowrbuf*2
            l=l+1
            r=spline%r_min*spline%p1**(i-1)
            if ( i .le. 0 ) then
                y(l)=spline%lb1*r + spline%lb2*r*r + spline%lb3*r*r*r + spline%lb4*r*r*r*r
            else
                f0=(i-1.0_r8)/real(lowrbuf*2,r8)
                f1=spline%lb1*r + spline%lb2*r*r + spline%lb3*r*r*r + spline%lb4*r*r*r*r
                y(l)=f0*coeff(1,i)+(1.0_r8-f0)*f1
            endif
            x(l)=i
        enddo
        ! Then tapered y-values
        do i=lowrbuf*2+1,n
            l=l+1
            ! store knot position
            x(l)=i
            ! Store value, with taper. Will update when I have more tapering thingies.
            if ( i .le. irc ) then
                y(l)=coeff(1,i)-yt(i,1)
            else
                y(l)=0.0_r8
            endif
        enddo

        ! Now construct the spline equations
        wy(1)=3*(y(2)-y(1))
        wd(1)=2
        wu(1)=1
        wl(1)=1
        do i=2,m-1
            wy(i)=3*(y(i+1)-y(i-1))
            wd(i)=4
            wl(i)=1
            wu(i)=1
        enddo
        ! And the final point
        wy(m)= 3*(y(m)-y(m-1))
        wd(m)= 2
        ! Solve tridiagonal system
        call dgtsv(n,1,wl,wd,wu,wy,n,i)
        if ( i .ne. 0 ) then
            call rl_stop_gracefully(['dgtsv exit code '//tochar(i)],rl_exitcode_blaslapack)
        endif
        ! Back-substitute spline parameters
        cleancoeff=0.0_r8
        do i=1,m-1
            cleancoeff(1,i) = y(i)
            cleancoeff(2,i) = wy(i)
            cleancoeff(3,i) = 3*(y(i+1)-y(i)) - 2*wy(i) - wy(i+1)
            cleancoeff(4,i) = 2*(y(i)-y(i+1)) + wy(i) + wy(i+1)
        enddo

        ! Copy only the necessary spline coefficients.
        spline%n_knot=0
        do i=1,size(coeff,2)
            r=spline%r_min*spline%p1**(i-1)
            if ( r .gt. cutoff+1E-5_r8 ) then
                exit
            else
                spline%n_knot=spline%n_knot+1
            endif
        enddo
        allocate(spline%coeff(4,spline%n_knot))
        spline%coeff=cleancoeff(:,lowrbuf+1:spline%n_knot+lowrbuf)
    !end block taperspline
    else
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
        allocate(spline%coeff_deriv(4,spline%n_knot))
        spline%coeff=coeff(:,1:spline%n_knot)
        spline%coeff_deriv=coeff_deriv(:,1:spline%n_knot)
    !end block donottaper
    endif
end subroutine

!> get the realspace rotation matrix for all basis functions on one atom
subroutine realspace_rotation_matrix_for_one_atom(basis,s1,op,rotationmatrix)
    !> basis set information
    class(rl_lcao_basis_set), intent(in) :: basis
    !> which species do we care about
    integer, intent(in) :: s1
    !> symmetry operation
    real(r8), dimension(3,3), intent(in) :: op
    !> rotation matrix
    real(r8), dimension(basis%species(s1)%n_basis,basis%species(s1)%n_basis), intent(out) :: rotationmatrix

    real(r8), dimension( 0:0, 0:0) :: rm0
    real(r8), dimension(-1:1,-1:1) :: rm1
    real(r8), dimension(-2:2,-2:2) :: rm2
    real(r8), dimension(-3:3,-3:3) :: rm3
    real(r8), dimension(-4:4,-4:4) :: rm4
    real(r8), dimension(-5:5,-5:5) :: rm5
    integer :: i1,i2,i,l

    ! If I use the inverse operation here, it actually means forward operation.
    rotationmatrix=0.0_r8
    do i=1,basis%species(s1)%n_basis_group
        l =basis%species(s1)%basis_group_l(i)
        i1=basis%species(s1)%basis_group(1,i)
        i2=basis%species(s1)%basis_group(2,i)
        select case(l)
        case(0)
            call olles_rotmat(op,0,rm0)
            rotationmatrix(i1:i2,i1:i2)=rm0(0:0,0:0)
        case(1)
            call olles_rotmat(op,1,rm1)
            rotationmatrix(i1:i2,i1:i2)=rm1(-1:1,-1:1)
        case(2)
            call olles_rotmat(op,2,rm2)
            rotationmatrix(i1:i2,i1:i2)=rm2(-2:2,-2:2)
        case(3)
            call olles_rotmat(op,3,rm3)
            rotationmatrix(i1:i2,i1:i2)=rm3(-3:3,-3:3)
        case(4)
            call olles_rotmat(op,4,rm4)
            rotationmatrix(i1:i2,i1:i2)=rm4(-4:4,-4:4)
        case(5)
            call olles_rotmat(op,5,rm5)
            rotationmatrix(i1:i2,i1:i2)=rm5(-5:5,-5:5)
        case default
            call rl_stop_gracefully(['Current max l is 5, nag on me to fix higher'],rl_exitcode_param)
        end select
    enddo
end subroutine

!> get the full transformation matrix for the basis functions, not sparse or anything
subroutine kspace_rotation_matrix(basis,op,p,qpoint,rotationmatrix)
    !> basis set information
    class(rl_lcao_basis_set), intent(in) :: basis
    !> symmetry operation
    type(rl_spacegroup_operation), intent(in) :: op
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> kpoint (in fractional Cartesian coordinates)
    real(r8), dimension(3), intent(in) :: qpoint
    !> rotation matrix
    complex(r8), dimension(basis%n_basis,basis%n_basis), intent(out) :: rotationmatrix

    real(r8), dimension(3,3) :: m0
    real(r8), dimension( 0:0, 0:0) :: rm0
    real(r8), dimension(-1:1,-1:1) :: rm1
    real(r8), dimension(-2:2,-2:2) :: rm2
    real(r8), dimension(-3:3,-3:3) :: rm3
    real(r8), dimension(-4:4,-4:4) :: rm4
    real(r8), dimension(-5:5,-5:5) :: rm5
    real(r8), dimension(3) :: v0,v1
    real(r8) :: arg
    complex(r8) :: phasefactor
    integer :: a1,a2,s1,s2,i1,i2,j1,j2,i,l

    ! Indexing of large matrix should be (new,old)
    ! Such that new = (new,old)*old
    ! If that makes any sense.
    m0=op%im
    rotationmatrix=0.0_r8
    do a1=1,p%n_atom
        a2=op%fmap(a1)
        s1=p%species(a1)
        s2=p%species(a2)
        do i=1,basis%species(s1)%n_basis_group
            ! Sort out where we are
            l =basis%species(s1)%basis_group_l(i)
            j1=basis%species(s1)%basis_group(1,i)+basis%offset(a1)
            j2=basis%species(s1)%basis_group(2,i)+basis%offset(a1)
            i1=basis%species(s2)%basis_group(1,i)+basis%offset(a2)
            i2=basis%species(s2)%basis_group(2,i)+basis%offset(a2)

            ! Get the phase factor. This is the same as from Maradudin.
            ! This actually works, tested brute force on very accurate
            ! eigenvectors, as opposed to the thingy in the AIMS manual.
            v0=qpoint*rl_twopi
            v1=matmul(op%im,p%cartesian_coordinate(:,a2))-p%cartesian_coordinate(:,a1)
            arg=dot_product(v0,v1)
            phasefactor=cmplx(cos(arg),sin(arg),r8)

            select case(l)
            case(0)
                call olles_rotmat(m0,0,rm0)
                rotationmatrix(i1:i2,j1:j2)=rm0(0:0,0:0)*phasefactor
            case(1)
                call olles_rotmat(m0,1,rm1)
                rotationmatrix(i1:i2,j1:j2)=rm1(-1:1,-1:1)*phasefactor
            case(2)
                call olles_rotmat(m0,2,rm2)
                rotationmatrix(i1:i2,j1:j2)=rm2(-2:2,-2:2)*phasefactor
            case(3)
                call olles_rotmat(m0,3,rm3)
                rotationmatrix(i1:i2,j1:j2)=rm3(-3:3,-3:3)*phasefactor
            case(4)
                call olles_rotmat(m0,4,rm4)
                rotationmatrix(i1:i2,j1:j2)=rm4(-4:4,-4:4)*phasefactor
            case(5)
                call olles_rotmat(m0,5,rm5)
                rotationmatrix(i1:i2,j1:j2)=rm5(-5:5,-5:5)*phasefactor
            case default
                call rl_stop_gracefully(['Current max l is 5, nag on me to fix higher'],rl_exitcode_param)
            end select
        enddo
    enddo
end subroutine

!> stupid rotation matrix I hate it but I think it works.
!@TODO store inversion and angles already in the symmetry operation
subroutine olles_rotmat(m,l,rotm)
    !> normal 3x3 rotation matrix
    real(r8), dimension(3,3), intent(in) :: m
    !> l number
    integer, intent(in) :: l
    !> rotation matrix
    real(r8), dimension(-l:l,-l:l), intent(out) :: rotm

    real(r8) :: al,be,gm
    real(r8), dimension(3,3) :: mm
    logical :: inversion

    !getrm: block
    real(r8) :: sa,sb,sg,ca,cb,cg,s2a,s2b,s2g,c2a,c2b,c2g,sa2,sb2,sg2,ca2,cb2,cg2
    real(r8) :: ca2sq,cb2sq,cg2sq,sa2sq,sb2sq,sg2sq
    real(r8) :: capg,camg,cbpg,cbmg,capb,camb,sapg,samg,sbpg,sbmg,sapb,samb
    real(r8), parameter :: half=0.5_r8,quarter=0.25_r8 !,threequarter=0.75_r8
    real(r8), parameter :: eigth=0.125_r8,inv16=1.0_r8/16.0_r8,inv32=1.0_r8/32.0_r8
    real(r8), parameter :: sqrt3=1.7320508075688772935_r8,sqrt5=2.2360679774997896964_r8,sqrt7=2.6457513110645905905_r8
    real(r8), parameter :: sqrt15=3.872983346207416885_r8,sqrt5o2=1.581138830084189666_r8,sqrt3o2=1.224744871391589049_r8
    real(r8), parameter :: sqrt35=5.916079783099616042_r8,invsqrt32=0.17677669529663688110_r8,inv64=0.015625_r8
    real(r8), parameter :: sqrt35o2=4.1833001326703777399_r8,sqrt7o2=1.8708286933869706928_r8,inv128=0.0078125_r8
    real(r8), parameter :: invsqrt128=0.088388347648318440550_r8,invsqrt2048=1.0_r8/sqrt(2048.0_r8)
    real(r8), parameter :: invsqrt512=1.0_r8/sqrt(512.0_r8),inv256=1.0_r8/256.0_r8,invsqrt2=1.0_r8/sqrt(2.0_r8)
    real(r8), parameter :: sqrt105=sqrt(105.0_r8),sqrt105o2=sqrt(105.0_r8/2.0_r8),sqrt15o2=sqrt(15.0_r8/2.0_r8)
    real(r8), parameter :: sqrt21=sqrt(21.0_r8),sqrt21o2=sqrt(21.0_r8/2.0_r8)

    ! check wether it's a proper or improper rotation
    if ( rl_determ(m) .lt. 0.0_r8 ) then
        mm=-m
        inversion=.true.
    else
        mm=m
        inversion=.false.
    endif

    ! grab the Euler angles for the proper rotation.
    ! this is some sort of convention. Not sure which.
    ! pretty sure this can be done in a more elegant way.
    if ( ( abs(mm(3,1)) .gt. rl_sqtol ) .or. ( abs(mm(3,2)) .gt. rl_sqtol ) ) then
        al = atan2(mm(3,2), mm(3,1))
        if ( al .lt. 0.0_r8 ) al=al+rl_twopi
        if ( abs(mm(3,1)) .gt. abs(mm(3,2)) ) then
            be = atan2(mm(3,1)/cos(al),mm(3,3))
        else
            be = atan2(mm(3,2)/sin(al),mm(3,3))
        endif
        gm = atan2(mm(2,3),-mm(1,3))
        if ( gm .lt. 0.0_r8 ) gm=gm+rl_twopi
    else
        al = atan2(mm(1,2),mm(1,1))
        if ( al .lt. 0.0_r8 ) al=al+rl_twopi
        if (mm(3,3) > 0d0) then
            be=0.0_r8
            gm=0.0_r8
        else
            be=rl_pi
            gm=rl_pi
        endif
    endif

    ! different for different l numbers. All tabulated neatly, up to l=5 for now. Not sure how high AIMS
    ! can go. Can fix more if needed, it's just a little annoying.
    !getrm: block
        ! a bunch of parameters and shorthand to make it a little faster.
        sa=sin(al)
        sb=sin(be)
        sg=sin(gm)
        ca=cos(al)
        cb=cos(be)
        cg=cos(gm)
        s2a=sin(2.0_r8*al)
        s2b=sin(2.0_r8*be)
        s2g=sin(2.0_r8*gm)
        c2a=cos(2.0_r8*al)
        c2b=cos(2.0_r8*be)
        c2g=cos(2.0_r8*gm)
        sa2=sin(0.5_r8*al)
        sb2=sin(0.5_r8*be)
        sg2=sin(0.5_r8*gm)
        ca2=cos(0.5_r8*al)
        cb2=cos(0.5_r8*be)
        cg2=cos(0.5_r8*gm)
        ca2sq=cos(al*0.5_r8)**2
        cb2sq=cos(be*0.5_r8)**2
        cg2sq=cos(gm*0.5_r8)**2
        sa2sq=sin(al*0.5_r8)**2
        sb2sq=sin(be*0.5_r8)**2
        sg2sq=sin(gm*0.5_r8)**2
        capg=cos(al+gm)
        camg=cos(al-gm)
        cbpg=cos(be+gm)
        cbmg=cos(be-gm)
        capb=cos(al+be)
        camb=cos(al-be)
        sapg=sin(al+gm)
        samg=sin(al-gm)
        sbpg=sin(be+gm)
        sbmg=sin(be-gm)
        sapb=sin(al+be)
        samb=sin(al-be)

        rotm=rl_huge
        select case(l)
        case(0) ! simple enough
            rotm=1.0_r8
        case(1) ! 3x3 thingy
            rotm(-1,-1)=capg*cb2sq + camg*sb2sq
            rotm( 0,-1)=sb*sg
            rotm( 1,-1)=cb2sq*sapg + samg*sb2sq
            rotm(-1, 0)=sa*sb
            rotm( 0, 0)=cb
            rotm( 1, 0)=-(ca*sb)
            rotm(-1, 1)=-(cb2sq*sapg) + samg*sb2sq
            rotm( 0, 1)=cg*sb
            rotm( 1, 1)=capg*cb2sq - camg*sb2sq
        case(2) ! 5x5 thingy
            rotm(-2,-2)=c2a*c2g*cb - (3 + c2b)*quarter*s2a*s2g
            rotm(-1,-2)=(-(c2g*ca) + cb*s2g*sa)*sb
            rotm(0,-2)=-(cg*sb**2*sg*sqrt3)
            rotm(1,-2)=-((ca*cb*s2g + c2g*sa)*sb)
            rotm(2,-2)=-(c2g*cb*s2a) - c2a*(3 + c2b)*quarter*s2g
            rotm(-2,-1)=sb*(c2a*cg - 2*ca*cb*sa*sg)
            rotm(-1,-1)=ca*cb*cg - c2b*sa*sg
            rotm(0,-1)=cb*sb*sg*sqrt3
            rotm(1,-1)=cb*cg*sa + c2b*ca*sg
            rotm(2,-1)=-(cg*s2a*sb) - c2a*half*s2b*sg
            rotm(-2,0)=ca*sa*sb**2*sqrt3
            rotm(-1,0)=cb*sa*sb*sqrt3
            rotm(0,0)=(1 + 3*c2b)*quarter
            rotm(1,0)=-(ca*cb*sb*sqrt3)
            rotm(2,0)=c2a*half*sb**2*sqrt3
            rotm(-2,1)=-(sb*(cb*cg*s2a + c2a*sg))
            rotm(-1,1)=-(c2b*cg*sa) - ca*cb*sg
            rotm(0,1)=cb*cg*sb*sqrt3
            rotm(1,1)=c2b*ca*cg - cb*sa*sg
            rotm(2,1)=-(c2a*cg*half*s2b) + s2a*sb*sg
            rotm(-2,2)=(3 + c2b)*c2g*quarter*s2a + c2a*cb*s2g
            rotm(-1,2)=-((ca*s2g + c2g*cb*sa)*sb)
            rotm(0,2)=c2g*half*sb**2*sqrt3
            rotm(1,2)=sb*(c2g*ca*cb - 2*cg*sa*sg)
            rotm(2,2)=c2a*(3 + c2b)*c2g*quarter - cb*s2a*s2g
        case(3) ! 7x7, getting hairy
            rotm(-3,-3)=sb2**6*Cos(3*(al - gm)) + cb2**6*Cos(3*(al + gm))
            rotm(-2,-3)=sb2*sqrt3o2*(sb*sb2**3*Cos(2*al - 3*gm) - 2*cb2**5*Cos(2*al + 3*gm))
            rotm(-1,-3)=quarter*sb2**2*sqrt15*(sb**2*Cos(al - 3*gm) + 4*cb2**4*Cos(al + 3*gm))
            rotm(0,-3)=half*sb**3*sqrt5o2*Sin(3*gm)
            rotm(1,-3)=quarter*sb2**2*sqrt15*(sb**2*Sin(al - 3*gm) + 4*cb2**4*Sin(al + 3*gm))
            rotm(2,-3)=sb2*sqrt3o2*(-(sb*sb2**3*Sin(2*al - 3*gm)) + 2*cb2**5*Sin(2*al + 3*gm))
            rotm(3,-3)=sb2**6*Sin(3*(al - gm)) + cb2**6*Sin(3*(al + gm))
            rotm(-3,-2)=sb2*sqrt3o2*(-(sb*sb2**3*Cos(3*al - 2*gm)) + 2*cb2**5*Cos(3*al + 2*gm))
            rotm(-2,-2)=-((2 + 3*cb)*sb2**4*Cos(2*(al - gm))) + (-2 + 3*cb)*cb2**4*Cos(2*(al + gm))
            rotm(-1,-2)=-(cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Cos(al - 2*gm) + (-1 + 3*cb)*cb2**2*Cos(al + 2*gm)))
            rotm(0,-2)=-(cb*cg*sb**2*sg*sqrt15)
            rotm(1,-2)=-(cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Sin(al - 2*gm) + (-1 + 3*cb)*cb2**2*Sin(al + 2*gm)))
            rotm(2,-2)=(2 + 3*cb)*sb2**4*Sin(2*(al - gm)) + (2 - 3*cb)*cb2**4*Sin(2*(al + gm))
            rotm(3,-2)=sb2*sqrt3o2*(-(sb*sb2**3*Sin(3*al - 2*gm)) + 2*cb2**5*Sin(3*al + 2*gm))
            rotm(-3,-1)=quarter*sb2**2*sqrt15*(sb**2*Cos(3*al - gm) + 4*cb2**4*Cos(3*al + gm))
            rotm(-2,-1)=cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Cos(2*al - gm) + (-1 + 3*cb)*cb2**2*Cos(2*al + gm))
            rotm(-1,-1)=eigth*((3 + 5*c2b)*ca*cg + (7 - 15*c2b)*cb*sa*sg)
            rotm(0,-1)=eigth*sg*sqrt3o2*(sb + 5*Sin(3*be))
            rotm(1,-1)=eigth*((3 + 5*c2b)*cg*sa + (-7 + 15*c2b)*ca*cb*sg)
            rotm(2,-1)=-(cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Sin(2*al - gm) + (-1 + 3*cb)*cb2**2*Sin(2*al + gm)))
            rotm(3,-1)=quarter*sb2**2*sqrt15*(sb**2*Sin(3*al - gm) + 4*cb2**4*Sin(3*al + gm))
            rotm(-3,0)=half*sb**3*sqrt5o2*Sin(3*al)
            rotm(-2,0)=ca*cb*sa*sb**2*sqrt15
            rotm(-1,0)=eigth*sa*sqrt3o2*(sb + 5*Sin(3*be))
            rotm(0,0)=eigth*(3*cb + 5*Cos(3*be))
            rotm(1,0)=-(ca*eigth*sqrt3o2*(sb + 5*Sin(3*be)))
            rotm(2,0)=c2a*cb*half*sb**2*sqrt15
            rotm(3,0)=-(half*sb**3*sqrt5o2*Cos(3*al))
            rotm(-3,1)=quarter*sb2**2*sqrt15*(sb**2*Sin(3*al - gm) - 4*cb2**4*Sin(3*al + gm))
            rotm(-2,1)=cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Sin(2*al - gm) + (1 - 3*cb)*cb2**2*Sin(2*al + gm))
            rotm(-1,1)=inv32*((43 - 45*c2b)*cb*cg*sa - 30*cb**3*cg*sa - 2*(1 + 5*c2b)*ca*sg - 20*ca*cb**2*sg)
            rotm(0,1)=cg*eigth*sqrt3o2*(sb + 5*Sin(3*be))
            rotm(1,1)=inv16*(-2*(3 + 5*c2b)*sa*sg + ca*cg*(cb + 15*Cos(3*be)))
            rotm(2,1)=cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Cos(2*al - gm) + (1 - 3*cb)*cb2**2*Cos(2*al + gm))
            rotm(3,1)=quarter*sb2**2*sqrt15*(-(sb**2*Cos(3*al - gm)) + 4*cb2**4*Cos(3*al + gm))
            rotm(-3,2)=sb2*sqrt3o2*(sb*sb2**3*Sin(3*al - 2*gm) + 2*cb2**5*Sin(3*al + 2*gm))
            rotm(-2,2)=(2 + 3*cb)*sb2**4*Sin(2*(al - gm)) + (-2 + 3*cb)*cb2**4*Sin(2*(al + gm))
            rotm(-1,2)=cb2*sb2*sqrt5o2*((1 + 3*cb)*sb2**2*Sin(al - 2*gm) + (1 - 3*cb)*cb2**2*Sin(al + 2*gm))
            rotm(0,2)=c2g*cb*half*sb**2*sqrt15
            rotm(1,2)=cb2*sb2*sqrt5o2*(-((1 + 3*cb)*sb2**2*Cos(al - 2*gm)) + (-1 + 3*cb)*cb2**2*Cos(al + 2*gm))
            rotm(2,2)=(2 + 3*cb)*sb2**4*Cos(2*(al - gm)) + (-2 + 3*cb)*cb2**4*Cos(2*(al + gm))
            rotm(3,2)=-(sb2*sqrt3o2*(sb*sb2**3*Cos(3*al - 2*gm) + 2*cb2**5*Cos(3*al + 2*gm)))
            rotm(-3,3)=sb2**6*Sin(3*(al - gm)) - cb2**6*Sin(3*(al + gm))
            rotm(-2,3)=sb2*sqrt3o2*(sb*sb2**3*Sin(2*al - 3*gm) + 2*cb2**5*Sin(2*al + 3*gm))
            rotm(-1,3)=quarter*sb2**2*sqrt15*(sb**2*Sin(al - 3*gm) - 4*cb2**4*Sin(al + 3*gm))
            rotm(0,3)=half*sb**3*sqrt5o2*Cos(3*gm)
            rotm(1,3)=quarter*sb2**2*sqrt15*(-(sb**2*Cos(al - 3*gm)) + 4*cb2**4*Cos(al + 3*gm))
            rotm(2,3)=sb2*sqrt3o2*(sb*sb2**3*Cos(2*al - 3*gm) + 2*cb2**5*Cos(2*al + 3*gm))
            rotm(3,3)=-(sb2**6*Cos(3*(al - gm))) + cb2**6*Cos(3*(al + gm))
        case(4) ! 9x9, getting rough.
            rotm(-4,-4)=eigth*Cos(4*al)*(7*cb + Cos(3*be))*Cos(4*gm) - inv64*(35 + 28*c2b + Cos(4*be))*Sin(4*al)*Sin(4*gm)
            rotm(-3,-4)=invsqrt512*(-2*Cos(3*al)*Cos(4*gm)*(7*sb + 3*Sin(3*be)) + Sin(3*al)*(14*s2b + Sin(4*be))*Sin(4*gm))
            rotm(-2,-4)=eigth*sb**2*sqrt7*(4*c2a*cb*Cos(4*gm) - (3 + c2b)*s2a*Sin(4*gm))
            rotm(-1,-4)=half*sb**3*sqrt7o2*(-(ca*Cos(4*gm)) + cb*sa*Sin(4*gm))
            rotm(0,-4)=-(eigth*sb**4*sqrt35*Sin(4*gm))
            rotm(1,-4)=-(half*sb**3*sqrt7o2*(sa*Cos(4*gm) + ca*cb*Sin(4*gm)))
            rotm(2,-4)=-(eigth*sb**2*sqrt7*(4*cb*s2a*Cos(4*gm) + c2a*(3 + c2b)*Sin(4*gm)))
            rotm(3,-4)=invsqrt2048*(-4*Cos(4*gm)*Sin(3*al)*(7*sb + 3*Sin(3*be)) - 2*Cos(3*al)*(14*s2b + Sin(4*be))*Sin(4*gm))
            rotm(4,-4)=-(eigth*(7*cb + Cos(3*be))*Cos(4*gm)*Sin(4*al)) - inv64*Cos(4*al)*(35 + 28*c2b + Cos(4*be))*Sin(4*gm)
            rotm(-4,-3)=invsqrt2048*(4*Cos(4*al)*Cos(3*gm)*(7*sb + 3*Sin(3*be)) - 2*Sin(4*al)*(14*s2b + Sin(4*be))*Sin(3*gm))
            rotm(-3,-3)=inv16*(Cos(3*al)*(7*cb + 9*Cos(3*be))*Cos(3*gm) - 2*(7*c2b + Cos(4*be))*Sin(3*al)*Sin(3*gm))
            rotm(-2,-3)=eigth*sqrt7o2*(c2a*Cos(3*gm)*(sb - 3*Sin(3*be)) + 8*cb**3*s2a*sb*Sin(3*gm))
            rotm(-1,-3)=quarter*sb**2*sqrt7*(3*ca*cb*Cos(3*gm) - (1 + 2*c2b)*sa*Sin(3*gm))
            rotm(0,-3)=cb*half*sb**3*sqrt35o2*Sin(3*gm)
            rotm(1,-3)=quarter*sb**2*sqrt7*(3*cb*sa*Cos(3*gm) + (1 + 2*c2b)*ca*Sin(3*gm))
            rotm(2,-3)=eigth*sqrt7o2*(s2a*Cos(3*gm)*(-sb + 3*Sin(3*be)) + 8*c2a*cb**3*sb*Sin(3*gm))
            rotm(3,-3)=inv16*((7*cb + 9*Cos(3*be))*Cos(3*gm)*Sin(3*al) + 2*Cos(3*al)*(7*c2b + Cos(4*be))*Sin(3*gm))
            rotm(4,-3)=invsqrt2048*(-4*Cos(3*gm)*Sin(4*al)*(7*sb + 3*Sin(3*be)) - 2*Cos(4*al)*(14*s2b + Sin(4*be))*Sin(3*gm))
            rotm(-4,-2)=eigth*sb**2*sqrt7*(4*c2g*cb*Cos(4*al) - (3 + c2b)*s2g*Sin(4*al))
            rotm(-3,-2)=-(eigth*sqrt7o2*(8*cb**3*s2g*sb*Sin(3*al) + c2g*Cos(3*al)*(sb - 3*Sin(3*be))))
            rotm(-2,-2)=c2a*c2g*eigth*(cb + 7*Cos(3*be)) - inv16*s2a*s2g*(5 + 4*c2b + 7*Cos(4*be))
            rotm(-1,-2)=-(invsqrt32*((5 + 7*c2b)*c2g*ca + 2*(1 - 7*c2b)*cb*s2g*sa)*sb)
            rotm(0,-2)=-((5 + 7*c2b)*eigth*s2g*sb**2*sqrt5)
            rotm(1,-2)=-(invsqrt128*(-2*ca*s2b*s2g + 3*c2g*sa*sb + 7*c2g*sa*Sin(3*be) + 7*ca*s2g*Sin(4*be)))
            rotm(2,-2)=-(c2g*eigth*s2a*(cb + 7*Cos(3*be))) - c2a*inv16*s2g*(5 + 4*c2b + 7*Cos(4*be))
            rotm(3,-2)=eigth*sqrt7o2*(8*cb**3*s2g*sb*Cos(3*al) + c2g*Sin(3*al)*(-sb + 3*Sin(3*be)))
            rotm(4,-2)=-(eigth*sb**2*sqrt7*((3 + c2b)*s2g*Cos(4*al) + 4*c2g*cb*Sin(4*al)))
            rotm(-4,-1)=half*sb**3*sqrt7o2*(cg*Cos(4*al) - cb*sg*Sin(4*al))
            rotm(-3,-1)=quarter*sb**2*sqrt7*(3*cb*cg*Cos(3*al) - (1 + 2*c2b)*sg*Sin(3*al))
            rotm(-2,-1)=invsqrt128*(3*c2a*cg*sb + 2*s2a*s2b*sg + 7*c2a*cg*Sin(3*be) - 7*s2a*sg*Sin(4*be))
            rotm(-1,-1)=inv16*(ca*cg*(9*cb + 7*Cos(3*be)) - 2*sa*sg*(c2b + 7*Cos(4*be)))
            rotm(0,-1)=(1 + 7*c2b)*eigth*s2b*sg*sqrt5o2
            rotm(1,-1)=inv16*(cg*sa*(9*cb + 7*Cos(3*be)) + 2*ca*sg*(c2b + 7*Cos(4*be)))
            rotm(2,-1)=-(invsqrt32*sb*((5 + 7*c2b)*cg*s2a + 5*c2a*cb*sg + 7*c2a*sg*Cos(3*be)))
            rotm(3,-1)=quarter*sb**2*sqrt7*((1 + 2*c2b)*sg*Cos(3*al) + 3*cb*cg*Sin(3*al))
            rotm(4,-1)=-(half*sb**3*sqrt7o2*(cb*sg*Cos(4*al) + cg*Sin(4*al)))
            rotm(-4,0)=eigth*sb**4*sqrt35*Sin(4*al)
            rotm(-3,0)=cb*half*sb**3*sqrt35o2*Sin(3*al)
            rotm(-2,0)=(5 + 7*c2b)*eigth*s2a*sb**2*sqrt5
            rotm(-1,0)=inv16*sa*sqrt5o2*(2*s2b + 7*Sin(4*be))
            rotm(0,0)=inv64*(9 + 20*c2b + 35*Cos(4*be))
            rotm(1,0)=-(ca*inv16*sqrt5o2*(2*s2b + 7*Sin(4*be)))
            rotm(2,0)=c2a*(5 + 7*c2b)*eigth*sb**2*sqrt5
            rotm(3,0)=-(cb*half*sb**3*sqrt35o2*Cos(3*al))
            rotm(4,0)=eigth*sb**4*sqrt35*Cos(4*al)
            rotm(-4,1)=-(half*sb**3*sqrt7o2*(sg*Cos(4*al) + cb*cg*Sin(4*al)))
            rotm(-3,1)=-(quarter*sb**2*sqrt7*(3*cb*sg*Cos(3*al) + (1 + 2*c2b)*cg*Sin(3*al)))
            rotm(-2,1)=-(invsqrt128*(-2*cg*s2a*s2b + 3*c2a*sb*sg + 7*c2a*sg*Sin(3*be) + 7*cg*s2a*Sin(4*be)))
            rotm(-1,1)=inv32*(-2*ca*sg*(9*cb + 7*Cos(3*be)) - 4*cg*sa*(c2b + 7*Cos(4*be)))
            rotm(0,1)=cg*inv16*sqrt5o2*(2*s2b + 7*Sin(4*be))
            rotm(1,1)=eigth*(-((1 + 7*c2b)*cb*sa*sg) + ca*cg*(c2b + 7*Cos(4*be)))
            rotm(2,1)=invsqrt128*(2*c2a*cg*s2b + 3*s2a*sb*sg + 7*s2a*sg*Sin(3*be) - 7*c2a*cg*Sin(4*be))
            rotm(3,1)=quarter*sb**2*sqrt7*((1 + 2*c2b)*cg*Cos(3*al) - 3*cb*sg*Sin(3*al))
            rotm(4,1)=half*sb**3*sqrt7o2*(-(cb*cg*Cos(4*al)) + sg*Sin(4*al))
            rotm(-4,2)=eigth*sb**2*sqrt7*(4*cb*s2g*Cos(4*al) + (3 + c2b)*c2g*Sin(4*al))
            rotm(-3,2)=eigth*sqrt7o2*(8*c2g*cb**3*sb*Sin(3*al) + s2g*Cos(3*al)*(-sb + 3*Sin(3*be)))
            rotm(-2,2)=inv16*(2*c2a*s2g*(cb + 7*Cos(3*be)) + c2g*s2a*(5 + 4*c2b + 7*Cos(4*be)))
            rotm(-1,2)=-(invsqrt128*(-2*c2g*s2b*sa + 3*ca*s2g*sb + 7*ca*s2g*Sin(3*be) + 7*c2g*sa*Sin(4*be)))
            rotm(0,2)=(5 + 7*c2b)*c2g*eigth*sb**2*sqrt5
            rotm(1,2)=invsqrt32*sb*(5*c2g*ca*cb - (5 + 7*c2b)*s2g*sa + 7*c2g*ca*Cos(3*be))
            rotm(2,2)=inv16*(-2*s2a*s2g*(cb + 7*Cos(3*be)) + c2a*c2g*(5 + 4*c2b + 7*Cos(4*be)))
            rotm(3,2)=-(eigth*sqrt7o2*(8*c2g*cb**3*sb*Cos(3*al) + s2g*Sin(3*al)*(sb - 3*Sin(3*be))))
            rotm(4,2)=eigth*sb**2*sqrt7*((3 + c2b)*c2g*Cos(4*al) - 4*cb*s2g*Sin(4*al))
            rotm(-4,3)=invsqrt2048*(-2*Cos(3*gm)*Sin(4*al)*(14*s2b + Sin(4*be)) - 4*Cos(4*al)*(7*sb + 3*Sin(3*be))*Sin(3*gm))
            rotm(-3,3)=-(eigth*(7*c2b + Cos(4*be))*Cos(3*gm)*Sin(3*al)) - inv16*Cos(3*al)*(7*cb + 9*Cos(3*be))*Sin(3*gm)
            rotm(-2,3)=eigth*sqrt7o2*(8*cb**3*s2a*sb*Cos(3*gm) + c2a*(-sb + 3*Sin(3*be))*Sin(3*gm))
            rotm(-1,3)=-(quarter*sb**2*sqrt7*((1 + 2*c2b)*sa*Cos(3*gm) + 3*ca*cb*Sin(3*gm)))
            rotm(0,3)=cb*half*sb**3*sqrt35o2*Cos(3*gm)
            rotm(1,3)=quarter*sb**2*sqrt7*((1 + 2*c2b)*ca*Cos(3*gm) - 3*cb*sa*Sin(3*gm))
            rotm(2,3)=eigth*sqrt7o2*(8*c2a*cb**3*sb*Cos(3*gm) + s2a*(sb - 3*Sin(3*be))*Sin(3*gm))
            rotm(3,3)=eigth*Cos(3*al)*(7*c2b + Cos(4*be))*Cos(3*gm) - inv16*(7*cb + 9*Cos(3*be))*Sin(3*al)*Sin(3*gm)
            rotm(4,3)=invsqrt2048*(-2*Cos(4*al)*Cos(3*gm)*(14*s2b + Sin(4*be)) + 4*Sin(4*al)*(7*sb + 3*Sin(3*be))*Sin(3*gm))
            rotm(-4,4)=inv64*((35 + 28*c2b + Cos(4*be))*Cos(4*gm)*Sin(4*al) + 8*Cos(4*al)*(7*cb + Cos(3*be))*Sin(4*gm))
            rotm(-3,4)=invsqrt2048*(-2*Cos(4*gm)*Sin(3*al)*(14*s2b + Sin(4*be)) - 4*Cos(3*al)*(7*sb + 3*Sin(3*be))*Sin(4*gm))
            rotm(-2,4)=eigth*sb**2*sqrt7*((3 + c2b)*s2a*Cos(4*gm) + 4*c2a*cb*Sin(4*gm))
            rotm(-1,4)=-(half*sb**3*sqrt7o2*(cb*sa*Cos(4*gm) + ca*Sin(4*gm)))
            rotm(0,4)=eigth*sb**4*sqrt35*Cos(4*gm)
            rotm(1,4)=half*sb**3*sqrt7o2*(ca*cb*Cos(4*gm) - sa*Sin(4*gm))
            rotm(2,4)=eigth*sb**2*sqrt7*(c2a*(3 + c2b)*Cos(4*gm) - 4*cb*s2a*Sin(4*gm))
            rotm(3,4)=invsqrt512*(Cos(3*al)*Cos(4*gm)*(14*s2b + Sin(4*be)) - 2*Sin(3*al)*(7*sb + 3*Sin(3*be))*Sin(4*gm))
            rotm(4,4)=inv64*(Cos(4*al)*(35 + 28*c2b + Cos(4*be))*Cos(4*gm) - 8*(7*cb + Cos(3*be))*Sin(4*al)*Sin(4*gm))
        case(5) ! 11x11, still works!
            rotm(-5,-5)=inv256*(2*Cos(5*al)*(63 + 60*c2b + 5*Cos(4*be))*Cos(5*gm) - (210*cb + 45*Cos(3*be) + Cos(5*be))*Sin(5*al)*Sin(5*gm))
            rotm(-4,-5)=inv128*sqrt5o2*(-8*Cos(4*al)*Cos(5*gm)*(6*s2b + Sin(4*be)) + Sin(4*al)*(42*sb + 27*Sin(3*be) + Sin(5*be))*Sin(5*gm))
            rotm(-3,-5)=3*inv64*sb**2*sqrt5*(2*(5 + 3*c2b)*Cos(3*al)*Cos(5*gm) - (15*cb + Cos(3*be))*Sin(3*al)*Sin(5*gm))
            rotm(-2,-5)=eigth*sb**3*sqrt15o2*(-4*c2a*cb*Cos(5*gm) + (3 + c2b)*s2a*Sin(5*gm))
            rotm(-1,-5)=eigth*sb**4*sqrt105o2*(ca*Cos(5*gm) - cb*sa*Sin(5*gm))
            rotm(0,-5)=3*eigth*sb**5*sg*sqrt7o2*(1 + 2*c2g + 2*Cos(4*gm))
            rotm(1,-5)=eigth*sb**4*sqrt105o2*(sa*Cos(5*gm) + ca*cb*Sin(5*gm))
            rotm(2,-5)=eigth*sb**3*sqrt15o2*(4*cb*s2a*Cos(5*gm) + c2a*(3 + c2b)*Sin(5*gm))
            rotm(3,-5)=3*inv64*sb**2*sqrt5*(2*(5 + 3*c2b)*Cos(5*gm)*Sin(3*al) + Cos(3*al)*(15*cb + Cos(3*be))*Sin(5*gm))
            rotm(4,-5)=inv128*sqrt5o2*(8*Cos(5*gm)*Sin(4*al)*(6*s2b + Sin(4*be)) + Cos(4*al)*(42*sb + 27*Sin(3*be) + Sin(5*be))*Sin(5*gm))
            rotm(5,-5)=inv256*(2*(63 + 60*c2b + 5*Cos(4*be))*Cos(5*gm)*Sin(5*al) + Cos(5*al)*(210*cb + 45*Cos(3*be) + Cos(5*be))*Sin(5*gm))
            rotm(-5,-4)=inv128*sqrt5o2*(8*Cos(5*al)*Cos(4*gm)*(6*s2b + Sin(4*be)) - Sin(5*al)*(42*sb + 27*Sin(3*be) + Sin(5*be))*Sin(4*gm))
            rotm(-4,-4)=quarter*Cos(4*al)*(3*c2b + Cos(4*be))*Cos(4*gm) - inv128*(42*cb + 81*Cos(3*be) + 5*Cos(5*be))*Sin(4*al)*Sin(4*gm)
            rotm(-3,-4)=-3*inv128*invsqrt2*(8*Cos(3*al)*Cos(4*gm)*(2*s2b + 3*Sin(4*be)) + Sin(3*al)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be))*Sin(4*gm))
            rotm(-2,-4)=inv16*sb**2*sqrt3*(8*c2a*(1 + 2*c2b)*Cos(4*gm) - s2a*(19*cb + 5*Cos(3*be))*Sin(4*gm))
            rotm(-1,-4)=inv16*sb**3*sqrt21*(-8*ca*cb*Cos(4*gm) + (3 + 5*c2b)*sa*Sin(4*gm))
            rotm(0,-4)=-3*cb*eigth*sb**4*sqrt35*Sin(4*gm)
            rotm(1,-4)=-(inv16*sb**3*sqrt21*(8*cb*sa*Cos(4*gm) + (3 + 5*c2b)*ca*Sin(4*gm)))
            rotm(2,-4)=-(inv16*sb**2*sqrt3*(8*(1 + 2*c2b)*s2a*Cos(4*gm) + c2a*(19*cb + 5*Cos(3*be))*Sin(4*gm)))
            rotm(3,-4)=-3*inv128*invsqrt2*(8*Cos(4*gm)*Sin(3*al)*(2*s2b + 3*Sin(4*be)) + Cos(3*al)*(-14*sb + 39*Sin(3*be) + 5*Sin(5*be))*Sin(4*gm))
            rotm(4,-4)=-(quarter*(3*c2b + Cos(4*be))*Cos(4*gm)*Sin(4*al)) - inv128*Cos(4*al)*(42*cb + 81*Cos(3*be) + 5*Cos(5*be))*Sin(4*gm)
            rotm(5,-4)=inv128*sqrt5o2*(8*Cos(4*gm)*Sin(5*al)*(6*s2b + Sin(4*be)) + Cos(5*al)*(42*sb + 27*Sin(3*be) + Sin(5*be))*Sin(4*gm))
            rotm(-5,-3)=3*inv64*sb**2*sqrt5*(2*(5 + 3*c2b)*Cos(5*al)*Cos(3*gm) - (15*cb + Cos(3*be))*Sin(5*al)*Sin(3*gm))
            rotm(-4,-3)=3*inv128*invsqrt2*(8*Cos(4*al)*Cos(3*gm)*(2*s2b + 3*Sin(4*be)) + Sin(4*al)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be))*Sin(3*gm))
            rotm(-3,-3)=inv256*(2*Cos(3*al)*(35 + 12*c2b + 81*Cos(4*be))*Cos(3*gm) - (42*cb + 169*Cos(3*be) + 45*Cos(5*be))*Sin(3*al)*Sin(3*gm))
            rotm(-2,-3)=inv64*sqrt3o2*(4*c2a*Cos(3*gm)*(2*s2b - 9*Sin(4*be)) + s2a*(14*sb + 13*Sin(3*be) + 15*Sin(5*be))*Sin(3*gm))
            rotm(-1,-3)=inv32*sb**2*sqrt21o2*(2*(7 + 9*c2b)*ca*Cos(3*gm) - sa*(17*cb + 15*Cos(3*be))*Sin(3*gm))
            rotm(0,-3)=(7 + 9*c2b)*inv16*sb**3*sqrt35o2*Sin(3*gm)
            rotm(1,-3)=inv32*sb**2*sqrt21o2*(2*(7 + 9*c2b)*sa*Cos(3*gm) + ca*(17*cb + 15*Cos(3*be))*Sin(3*gm))
            rotm(2,-3)=inv64*sqrt3o2*(4*s2a*Cos(3*gm)*(-2*s2b + 9*Sin(4*be)) + c2a*(14*sb + 13*Sin(3*be) + 15*Sin(5*be))*Sin(3*gm))
            rotm(3,-3)=inv256*(2*(35 + 12*c2b + 81*Cos(4*be))*Cos(3*gm)*Sin(3*al) + Cos(3*al)*(42*cb + 169*Cos(3*be) + 45*Cos(5*be))*Sin(3*gm))
            rotm(4,-3)=-3*inv128*invsqrt2*(8*Cos(3*gm)*Sin(4*al)*(2*s2b + 3*Sin(4*be)) + Cos(4*al)*(-14*sb + 39*Sin(3*be) + 5*Sin(5*be))*Sin(3*gm))
            rotm(5,-3)=3*inv64*sb**2*sqrt5*(2*(5 + 3*c2b)*Cos(3*gm)*Sin(5*al) + Cos(5*al)*(15*cb + Cos(3*be))*Sin(3*gm))
            rotm(-5,-2)=eigth*sb**3*sqrt15o2*(4*c2g*cb*Cos(5*al) - (3 + c2b)*s2g*Sin(5*al))
            rotm(-4,-2)=inv16*sb**2*sqrt3*(8*(1 + 2*c2b)*c2g*Cos(4*al) - s2g*(19*cb + 5*Cos(3*be))*Sin(4*al))
            rotm(-3,-2)=-(inv64*sqrt3o2*(4*c2g*Cos(3*al)*(2*s2b - 9*Sin(4*be)) + s2g*Sin(3*al)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be))))
            rotm(-2,-2)=c2a*c2g*quarter*(c2b + 3*Cos(4*be)) - inv32*s2a*s2g*(14*cb + 3*Cos(3*be) + 15*Cos(5*be))
            rotm(-1,-2)=-(inv32*sb*sqrt7*(4*c2g*ca*(5*cb + 3*Cos(3*be)) - s2g*sa*(5 + 12*c2b + 15*Cos(4*be))))
            rotm(0,-2)=-((1 + 3*c2b)*cb*eigth*s2g*sb**2*sqrt105)
            rotm(1,-2)=inv64*sqrt7*(-4*c2g*sa*(2*s2b + 3*Sin(4*be)) + ca*s2g*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(2,-2)=-(c2g*quarter*s2a*(c2b + 3*Cos(4*be))) - c2a*inv32*s2g*(14*cb + 3*Cos(3*be) + 15*Cos(5*be))
            rotm(3,-2)=inv64*sqrt3o2*(4*c2g*Sin(3*al)*(-2*s2b + 9*Sin(4*be)) + s2g*Cos(3*al)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be)))
            rotm(4,-2)=-(inv16*sb**2*sqrt3*(s2g*Cos(4*al)*(19*cb + 5*Cos(3*be)) + 8*(1 + 2*c2b)*c2g*Sin(4*al)))
            rotm(5,-2)=eigth*sb**3*sqrt15o2*((3 + c2b)*s2g*Cos(5*al) + 4*c2g*cb*Sin(5*al))
            rotm(-5,-1)=eigth*sb**4*sqrt105o2*(cg*Cos(5*al) - cb*sg*Sin(5*al))
            rotm(-4,-1)=inv16*sb**3*sqrt21*(8*cb*cg*Cos(4*al) - (3 + 5*c2b)*sg*Sin(4*al))
            rotm(-3,-1)=inv32*sb**2*sqrt21o2*(2*(7 + 9*c2b)*cg*Cos(3*al) - sg*(17*cb + 15*Cos(3*be))*Sin(3*al))
            rotm(-2,-1)=inv64*sqrt7*(4*c2a*cg*(2*s2b + 3*Sin(4*be)) + s2a*sg*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(-1,-1)=inv64*(cb*sa*sg*(-43 + 84*c2b - 105*Cos(4*be)) + ca*cg*(15 + 28*c2b + 21*Cos(4*be)))
            rotm(0,-1)=inv64*sb*sg*sqrt15*(15 + 28*c2b + 21*Cos(4*be))
            rotm(1,-1)=inv64*(cg*sa*(15 + 28*c2b + 21*Cos(4*be)) + ca*cb*sg*(43 - 84*c2b + 105*Cos(4*be)))
            rotm(2,-1)=inv64*sqrt7*(-4*cg*s2a*(2*s2b + 3*Sin(4*be)) + c2a*sg*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(3,-1)=inv32*sb**2*sqrt21o2*(sg*Cos(3*al)*(17*cb + 15*Cos(3*be)) + 2*(7 + 9*c2b)*cg*Sin(3*al))
            rotm(4,-1)=-(inv16*sb**3*sqrt21*((3 + 5*c2b)*sg*Cos(4*al) + 8*cb*cg*Sin(4*al)))
            rotm(5,-1)=eigth*sb**4*sqrt105o2*(cb*sg*Cos(5*al) + cg*Sin(5*al))
            rotm(-5,0)=3*eigth*sa*sb**5*sqrt7o2*(1 + 2*c2a + 2*Cos(4*al))
            rotm(-4,0)=3*cb*eigth*sb**4*sqrt35*Sin(4*al)
            rotm(-3,0)=(1 + 2*c2a)*(7 + 9*c2b)*inv16*sa*sb**3*sqrt35o2
            rotm(-2,0)=(1 + 3*c2b)*cb*eigth*s2a*sb**2*sqrt105
            rotm(-1,0)=inv64*sa*sb*sqrt15*(15 + 28*c2b + 21*Cos(4*be))
            rotm(0,0)=inv128*(30*cb + 35*Cos(3*be) + 63*Cos(5*be))
            rotm(1,0)=-(ca*inv128*sqrt15*(2*sb + 7*Sin(3*be) + 21*Sin(5*be)))
            rotm(2,0)=c2a*inv16*sb**2*sqrt105*(5*cb + 3*Cos(3*be))
            rotm(3,0)=-((7 + 9*c2b)*inv16*sb**3*sqrt35o2*Cos(3*al))
            rotm(4,0)=3*cb*eigth*sb**4*sqrt35*Cos(4*al)
            rotm(5,0)=-3*eigth*sb**5*sqrt7o2*Cos(5*al)
            rotm(-5,1)=-(eigth*sb**4*sqrt105o2*(sg*Cos(5*al) + cb*cg*Sin(5*al)))
            rotm(-4,1)=-(inv16*sb**3*sqrt21*(8*cb*sg*Cos(4*al) + (3 + 5*c2b)*cg*Sin(4*al)))
            rotm(-3,1)=-(inv32*sb**2*sqrt21o2*(2*(7 + 9*c2b)*sg*Cos(3*al) + cg*(17*cb + 15*Cos(3*be))*Sin(3*al)))
            rotm(-2,1)=inv64*sqrt7*(-4*c2a*sg*(2*s2b + 3*Sin(4*be)) + cg*s2a*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(-1,1)=inv128*(-2*ca*sg*(15 + 28*c2b + 21*Cos(4*be)) - cg*sa*(2*cb + 21*(Cos(3*be) + 5*Cos(5*be))))
            rotm(0,1)=cg*inv128*sqrt15*(2*sb + 7*Sin(3*be) + 21*Sin(5*be))
            rotm(1,1)=inv128*(-2*sa*sg*(15 + 28*c2b + 21*Cos(4*be)) + ca*cg*(2*cb + 21*(Cos(3*be) + 5*Cos(5*be))))
            rotm(2,1)=inv64*sqrt7*(4*s2a*sg*(2*s2b + 3*Sin(4*be)) + c2a*cg*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(3,1)=inv32*sb**2*sqrt21o2*(cg*Cos(3*al)*(17*cb + 15*Cos(3*be)) - 2*(7 + 9*c2b)*sg*Sin(3*al))
            rotm(4,1)=inv16*sb**3*sqrt21*(-((3 + 5*c2b)*cg*Cos(4*al)) + 8*cb*sg*Sin(4*al))
            rotm(5,1)=eigth*sb**4*sqrt105o2*(cb*cg*Cos(5*al) - sg*Sin(5*al))
            rotm(-5,2)=eigth*sb**3*sqrt15o2*(4*cb*s2g*Cos(5*al) + (3 + c2b)*c2g*Sin(5*al))
            rotm(-4,2)=inv16*sb**2*sqrt3*(8*(1 + 2*c2b)*s2g*Cos(4*al) + c2g*(19*cb + 5*Cos(3*be))*Sin(4*al))
            rotm(-3,2)=inv64*sqrt3o2*(4*s2g*Cos(3*al)*(-2*s2b + 9*Sin(4*be)) + c2g*Sin(3*al)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be)))
            rotm(-2,2)=inv32*(8*c2a*s2g*(c2b + 3*Cos(4*be)) + c2g*s2a*(14*cb + 3*Cos(3*be) + 15*Cos(5*be)))
            rotm(-1,2)=inv64*sqrt7*(-4*ca*s2g*(2*s2b + 3*Sin(4*be)) + c2g*sa*(2*sb + 3*Sin(3*be) - 15*Sin(5*be)))
            rotm(0,2)=(1 + 3*c2b)*c2g*cb*eigth*sb**2*sqrt105
            rotm(1,2)=inv32*sb*sqrt7*(-4*s2g*sa*(5*cb + 3*Cos(3*be)) + c2g*ca*(5 + 12*c2b + 15*Cos(4*be)))
            rotm(2,2)=inv32*(-8*s2a*s2g*(c2b + 3*Cos(4*be)) + c2a*c2g*(14*cb + 3*Cos(3*be) + 15*Cos(5*be)))
            rotm(3,2)=-(inv64*sqrt3o2*(4*s2g*Sin(3*al)*(2*s2b - 9*Sin(4*be)) + c2g*Cos(3*al)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be))))
            rotm(4,2)=inv16*sb**2*sqrt3*(c2g*Cos(4*al)*(19*cb + 5*Cos(3*be)) - 8*(1 + 2*c2b)*s2g*Sin(4*al))
            rotm(5,2)=eigth*sb**3*sqrt15o2*(-((3 + c2b)*c2g*Cos(5*al)) + 4*cb*s2g*Sin(5*al))
            rotm(-5,3)=-3*inv64*sb**2*sqrt5*((15*cb + Cos(3*be))*Cos(3*gm)*Sin(5*al) + 2*(5 + 3*c2b)*Cos(5*al)*Sin(3*gm))
            rotm(-4,3)=3*inv128*invsqrt2*(Cos(3*gm)*Sin(4*al)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be)) - 8*Cos(4*al)*(2*s2b + 3*Sin(4*be))*Sin(3*gm))
            rotm(-3,3)=inv256*(-((42*cb + 169*Cos(3*be) + 45*Cos(5*be))*Cos(3*gm)*Sin(3*al)) - 2*Cos(3*al)*(35 + 12*c2b + 81*Cos(4*be))*Sin(3*gm))
            rotm(-2,3)=inv64*sqrt3o2*(s2a*Cos(3*gm)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be)) + 4*c2a*(-2*s2b + 9*Sin(4*be))*Sin(3*gm))
            rotm(-1,3)=-(inv32*sb**2*sqrt21o2*(sa*(17*cb + 15*Cos(3*be))*Cos(3*gm) + 2*(7 + 9*c2b)*ca*Sin(3*gm)))
            rotm(0,3)=(7 + 9*c2b)*inv16*sb**3*sqrt35o2*Cos(3*gm)
            rotm(1,3)=inv32*sb**2*sqrt21o2*(ca*(17*cb + 15*Cos(3*be))*Cos(3*gm) - 2*(7 + 9*c2b)*sa*Sin(3*gm))
            rotm(2,3)=inv64*sqrt3o2*(c2a*Cos(3*gm)*(14*sb + 13*Sin(3*be) + 15*Sin(5*be)) + 4*s2a*(2*s2b - 9*Sin(4*be))*Sin(3*gm))
            rotm(3,3)=inv256*(Cos(3*al)*(42*cb + 169*Cos(3*be) + 45*Cos(5*be))*Cos(3*gm) - 2*(35 + 12*c2b + 81*Cos(4*be))*Sin(3*al)*Sin(3*gm))
            rotm(4,3)=3*inv128*invsqrt2*(Cos(4*al)*Cos(3*gm)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be)) + 8*Sin(4*al)*(2*s2b + 3*Sin(4*be))*Sin(3*gm))
            rotm(5,3)=3*inv64*sb**2*sqrt5*(Cos(5*al)*(15*cb + Cos(3*be))*Cos(3*gm) - 2*(5 + 3*c2b)*Sin(5*al)*Sin(3*gm))
            rotm(-5,4)=inv128*sqrt5o2*(Cos(4*gm)*Sin(5*al)*(42*sb + 27*Sin(3*be) + Sin(5*be)) + 8*Cos(5*al)*(6*s2b + Sin(4*be))*Sin(4*gm))
            rotm(-4,4)=inv128*(42*cb + 81*Cos(3*be) + 5*Cos(5*be))*Cos(4*gm)*Sin(4*al) + quarter*Cos(4*al)*(3*c2b + Cos(4*be))*Sin(4*gm)
            rotm(-3,4)=3*inv128*invsqrt2*(Cos(4*gm)*Sin(3*al)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be)) - 8*Cos(3*al)*(2*s2b + 3*Sin(4*be))*Sin(4*gm))
            rotm(-2,4)=inv16*sb**2*sqrt3*(s2a*(19*cb + 5*Cos(3*be))*Cos(4*gm) + 8*c2a*(1 + 2*c2b)*Sin(4*gm))
            rotm(-1,4)=-(inv16*sb**3*sqrt21*((3 + 5*c2b)*sa*Cos(4*gm) + 8*ca*cb*Sin(4*gm)))
            rotm(0,4)=3*cb*eigth*sb**4*sqrt35*Cos(4*gm)
            rotm(1,4)=inv16*sb**3*sqrt21*((3 + 5*c2b)*ca*Cos(4*gm) - 8*cb*sa*Sin(4*gm))
            rotm(2,4)=inv16*sb**2*sqrt3*(c2a*(19*cb + 5*Cos(3*be))*Cos(4*gm) - 8*(1 + 2*c2b)*s2a*Sin(4*gm))
            rotm(3,4)=-3*inv128*invsqrt2*(Cos(3*al)*Cos(4*gm)*(14*sb - 39*Sin(3*be) - 5*Sin(5*be)) + 8*Sin(3*al)*(2*s2b + 3*Sin(4*be))*Sin(4*gm))
            rotm(4,4)=inv128*Cos(4*al)*(42*cb + 81*Cos(3*be) + 5*Cos(5*be))*Cos(4*gm) - quarter*(3*c2b + Cos(4*be))*Sin(4*al)*Sin(4*gm)
            rotm(5,4)=inv128*sqrt5o2*(-(Cos(5*al)*Cos(4*gm)*(42*sb + 27*Sin(3*be) + Sin(5*be))) + 8*Sin(5*al)*(6*s2b + Sin(4*be))*Sin(4*gm))
            rotm(-5,5)=inv256*(-((210*cb + 45*Cos(3*be) + Cos(5*be))*Cos(5*gm)*Sin(5*al)) - 2*Cos(5*al)*(63 + 60*c2b + 5*Cos(4*be))*Sin(5*gm))
            rotm(-4,5)=inv128*sqrt5o2*(Cos(5*gm)*Sin(4*al)*(42*sb + 27*Sin(3*be) + Sin(5*be)) + 8*Cos(4*al)*(6*s2b + Sin(4*be))*Sin(5*gm))
            rotm(-3,5)=-3*inv64*sb**2*sqrt5*((15*cb + Cos(3*be))*Cos(5*gm)*Sin(3*al) + 2*(5 + 3*c2b)*Cos(3*al)*Sin(5*gm))
            rotm(-2,5)=eigth*sb**3*sqrt15o2*((3 + c2b)*s2a*Cos(5*gm) + 4*c2a*cb*Sin(5*gm))
            rotm(-1,5)=-(eigth*sb**4*sqrt105o2*(cb*sa*Cos(5*gm) + ca*Sin(5*gm)))
            rotm(0,5)=3*eigth*sb**5*sqrt7o2*Cos(5*gm)
            rotm(1,5)=eigth*sb**4*sqrt105o2*(ca*cb*Cos(5*gm) - sa*Sin(5*gm))
            rotm(2,5)=eigth*sb**3*sqrt15o2*(c2a*(3 + c2b)*Cos(5*gm) - 4*cb*s2a*Sin(5*gm))
            rotm(3,5)=3*inv64*sb**2*sqrt5*(Cos(3*al)*(15*cb + Cos(3*be))*Cos(5*gm) - 2*(5 + 3*c2b)*Sin(3*al)*Sin(5*gm))
            rotm(4,5)=inv128*sqrt5o2*(Cos(4*al)*Cos(5*gm)*(42*sb + 27*Sin(3*be) + Sin(5*be)) - 8*Sin(4*al)*(6*s2b + Sin(4*be))*Sin(5*gm))
            rotm(5,5)=inv256*(Cos(5*al)*(210*cb + 45*Cos(3*be) + Cos(5*be))*Cos(5*gm) - 2*(63 + 60*c2b + 5*Cos(4*be))*Sin(5*al)*Sin(5*gm))
        case default
        end select
    !end block getrm

    ! fix inversion
    if ( inversion ) then
        if ( mod(l,2) .ne. 0 ) rotm=-rotm
    endif
    ! and slice away tiny elements, set zeros to zero and ones to one.
    rotm=rl_chop(rotm,1E-11_r8)
end subroutine

!> evaluate radial part of basis function
elemental subroutine radial_function(h,r,y)
    !> function handle
    class(rl_lcao_basis_set_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y

    integer :: i
    real(r8) :: a0,a1,a2,a3,alpha

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        return
    elseif ( r .lt. h%r_min ) then
        y=h%lb1 + h%lb2*r + h%lb3*r*r + h%lb4*r*r*r
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
    ! See if sending in ir gives a speedup.
    y=y/r
end subroutine

!> evaluate radial part of basis function and the derivative
elemental subroutine radial_function_derivative(h,r,y,dy)
    !> function handle
    class(rl_lcao_basis_set_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y
    !> derivative of radial function
    real(r8), intent(out) :: dy

    integer :: i,j
    real(r8) :: a0,a1,a2,a3,alpha
    real(r8) :: i2,i3,yalpha,dydalpha,ir
    real(r8) :: x

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        dy=0.0_r8
        return
    elseif ( r .lt. h%r_min ) then
        y=h%lb1 + h%lb2*r + h%lb3*r*r + h%lb4*r*r*r
        dy=h%lb2 + 2*h%lb3*r + 3*h%lb4*r*r
        return
    endif

    ! Invert log grid coordinates
    !r = 1 + log(rad/h%smallest_radial_point)/log(h%radial_increment)
    alpha = 1 + log(r*h%ip0)*h%ilp1 ! Coordinate in index-space
    i = int(alpha)                  ! Index to spline
    alpha=alpha-i                   ! Shifted coordinate in index-space
    a0=h%coeff(1,i)
    a1=h%coeff(2,i)
    a2=h%coeff(3,i)
    a3=h%coeff(4,i)

    !i2=a2*alpha
    !i3=a3*alpha*alpha
    !
    !yalpha = a0 + a1*alpha + i2*alpha + i3*alpha
    !dydalpha = a1 + 2*i2 + 3*i3

    yalpha = a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha

    a0=h%coeff_deriv(1,i)
    a1=h%coeff_deriv(2,i)
    a2=h%coeff_deriv(3,i)
    a3=h%coeff_deriv(4,i)

    dydalpha = a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha

    ! Function value
    ir=1.0_r8/r
    y=yalpha*ir
    ! Derivative
    !dy=dydalpha !*ir !(dydalpha*h%ilp1-yalpha)*(ir*ir)
    dy=dydalpha*ir-yalpha*ir*ir
end subroutine

!> evaluate radial part of basis function and the derivative
elemental subroutine radial_function_derivative_hessian(h,r,y,dy,ddy)
    !> function handle
    class(rl_lcao_basis_set_spline), intent(in) :: h
    !> distance of point
    real(r8), intent(in) :: r
    !> function value
    real(r8), intent(out) :: y
    !> derivative of radial function
    real(r8), intent(out) :: dy
    !> second derivative of radial function
    real(r8), intent(out) :: ddy

    integer :: i
    real(r8) :: a0,a1,a2,a3,alpha
    real(r8) :: i2,i3,j3,yalpha,dydalpha,ddydalpha,ir

    ! Check if we are in any if the special regions.
    if ( r .ge. h%r_max ) then
        y=0.0_r8
        dy=0.0_r8
        ddy=0.0_r8
        return
    elseif ( r .lt. h%r_min ) then
        y=h%lb1 + h%lb2*r + h%lb3*r*r + h%lb4*r*r*r
        dy=h%lb2 + 2*h%lb3*r + 3*h%lb4*r*r
        ddy=2*h%lb3 + 6*h%lb4*r
        return
    endif

    ! Invert log grid coordinates
    !r = 1 + log(rad/h%smallest_radial_point)/log(h%radial_increment)
    alpha = 1 + log(r*h%ip0)*h%ilp1 ! Coordinate in index-space
    i = int(alpha)                  ! Index to spline
    alpha=alpha-i                   ! Shifted coordinate in index-space
    a0=h%coeff(1,i)
    a1=h%coeff(2,i)
    a2=h%coeff(3,i)
    a3=h%coeff(4,i)

    !  y  = a0 + a1*alpha +   a2*alpha**2 +   a3*alpha**3
    ! dy  =      a1       + 2*a2*alpha    + 3*a3*alpha**2
    ! ddy =                 2*a2          + 6*a3*alpha
    j3=a3*alpha
    i3=j3*alpha
    i2=a2*alpha

    ddydalpha = 2*a2+6*j3
    dydalpha = a1 + 2*i2+3*i3
    yalpha = a0 + a1*alpha + i2*alpha + i3*alpha

    ! Function value
    ir=1.0_r8/r
    y=yalpha*ir
    ! Derivative
    dy=(dydalpha*h%ilp1-yalpha)*(ir*ir)
    ! Second derivative
    ddy=(2*yalpha-3*dydalpha*h%ilp1+ddydalpha*h%ilp1*h%ilp1)*(ir*ir*ir)
end subroutine

!> evaluate spherical harmonic functions
elemental subroutine spherical_harmonic(xh,yh,zh,iang,y)
    !> xhat, xh=x/r
    real(r8), intent(in) :: xh
    !> yhat, yh=y/r
    real(r8), intent(in) :: yh
    !> zhat, zh=z/r
    real(r8), intent(in) :: zh
    !> angular index, iang=l^2+l+m+1
    integer, intent(in) :: iang
    !> spherical harmonic
    real(r8), intent(out) :: y

    ! [xh,yh,zh] is a unit vector, or in other words
    ! xh=x/r and so on. The l and m indices are contracted to
    ! a single index. For now, the spherical harmonics are only
    ! defined up to l=6, but it's easy to add more as needed.
    select case(iang)
    case(1)
        ! l=0, m=0, iang=1 (=l^2+l+1+m)
        ! Ylm= 1/(2*Sqrt[Pi])
        y=2.820947917738781E-1_r8
    case(2)
        ! l=1, m=-1, iang=2 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*yh)/2
        y=4.886025119029199E-1_r8*(yh)
    case(3)
        ! l=1, m=0, iang=3 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*zh)/2
        y=4.886025119029199E-1_r8*(zh)
    case(4)
        ! l=1, m=1, iang=4 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[3/Pi]*xh)/2
        y=-4.886025119029199E-1_r8*(xh)
    case(5)
        ! l=2, m=-2, iang=5 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*xh*yh)/2
        y=1.092548430592079E0_r8*(xh*yh)
    case(6)
        ! l=2, m=-1, iang=6 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*yh*zh)/2
        y=1.092548430592079E0_r8*(yh*zh)
    case(7)
        ! l=2, m=0, iang=7 (=l^2+l+1+m)
        ! Ylm= (Sqrt[5/Pi]*(-1 + 3*zh^2))/4
        y=3.153915652525200E-1_r8*(-1 + 3*zh**2)
    case(8)
        ! l=2, m=1, iang=8 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[15/Pi]*xh*zh)/2
        y=-1.092548430592079E0_r8*(xh*zh)
    case(9)
        ! l=2, m=2, iang=9 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*(xh^2 - yh^2))/4
        y=5.462742152960395E-1_r8*(xh**2 - yh**2)
    case(10)
        ! l=3, m=-3, iang=10 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2))/4
        y=-5.900435899266435E-1_r8*(-3*xh**2*yh + yh**3)
    case(11)
        ! l=3, m=-2, iang=11 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*xh*yh*zh)/2
        y=2.890611442640554E0_r8*(xh*yh*zh)
    case(12)
        ! l=3, m=-1, iang=12 (=l^2+l+1+m)
        ! Ylm= (Sqrt[21/(2*Pi)]*yh*(-1 + 5*zh^2))/4
        y=4.570457994644657E-1_r8*(yh*(-1 + 5*zh**2))
    case(13)
        ! l=3, m=0, iang=13 (=l^2+l+1+m)
        ! Ylm= (Sqrt[7/Pi]*zh*(-3 + 5*zh^2))/4
        y=3.731763325901154E-1_r8*(zh*(-3 + 5*zh**2))
    case(14)
        ! l=3, m=1, iang=14 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[21/(2*Pi)]*xh*(-1 + 5*zh^2))/4
        y=-4.570457994644657E-1_r8*(xh*(-1 + 5*zh**2))
    case(15)
        ! l=3, m=2, iang=15 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*(xh^2 - yh^2)*zh)/4
        y=1.445305721320277E0_r8*((xh**2 - yh**2)*zh)
    case(16)
        ! l=3, m=3, iang=16 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2))/4
        y=-5.900435899266435E-1_r8*(xh**3 - 3*xh*yh**2)
    case(17)
        ! l=4, m=-4, iang=17 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*xh*yh*(xh^2 - yh^2))/4
        y=2.503342941796705E0_r8*(xh*yh*(xh**2 - yh**2))
    case(18)
        ! l=4, m=-3, iang=18 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh)/4
        y=-1.770130769779931E0_r8*(yh*(-3*xh**2 + yh**2)*zh)
    case(19)
        ! l=4, m=-2, iang=19 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*xh*yh*(-1 + 7*zh^2))/4
        y=9.461746957575600E-1_r8*(xh*yh*(-1 + 7*zh**2))
    case(20)
        ! l=4, m=-1, iang=20 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/(2*Pi)]*yh*zh*(-3 + 7*zh^2))/4
        y=6.690465435572892E-1_r8*(yh*zh*(-3 + 7*zh**2))
    case(21)
        ! l=4, m=0, iang=21 (=l^2+l+1+m)
        ! Ylm= (3*(3 - 30*zh^2 + 35*zh^4))/(16*Sqrt[Pi])
        y=1.057855469152043E-1_r8*(3 - 30*zh**2 + 35*zh**4)
    case(22)
        ! l=4, m=1, iang=22 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[5/(2*Pi)]*xh*zh*(-3 + 7*zh^2))/4
        y=-6.690465435572892E-1_r8*(xh*zh*(-3 + 7*zh**2))
    case(23)
        ! l=4, m=2, iang=23 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*(xh^2 - yh^2)*(-1 + 7*zh^2))/8
        y=4.730873478787800E-1_r8*((xh**2 - yh**2)*(-1 + 7*zh**2))
    case(24)
        ! l=4, m=3, iang=24 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh)/4
        y=-1.770130769779931E0_r8*(xh*(xh**2 - 3*yh**2)*zh)
    case(25)
        ! l=4, m=4, iang=25 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4))/16
        y=6.258357354491761E-1_r8*(xh**4 - 6*xh**2*yh**2 + yh**4)
    case(26)
        ! l=5, m=-5, iang=26 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4))/16
        y=6.563820568401701E-1_r8*(5*xh**4*yh - 10*xh**2*yh**3 + yh**5)
    case(27)
        ! l=5, m=-4, iang=27 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*xh*yh*(xh^2 - yh^2)*zh)/4
        y=8.302649259524165E0_r8*(xh*yh*(xh**2 - yh**2)*zh)
    case(28)
        ! l=5, m=-3, iang=28 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 9*zh^2))/16
        y=-4.892382994352504E-1_r8*(yh*(-3*xh**2 + yh**2)*(-1 + 9*zh**2))
    case(29)
        ! l=5, m=-2, iang=29 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*xh*yh*zh*(-1 + 3*zh^2))/4
        y=4.793536784973324E0_r8*(xh*yh*zh*(-1 + 3*zh**2))
    case(30)
        ! l=5, m=-1, iang=30 (=l^2+l+1+m)
        ! Ylm= (Sqrt[165/Pi]*yh*(1 - 14*zh^2 + 21*zh^4))/16
        y=4.529466511956969E-1_r8*(yh*(1 - 14*zh**2 + 21*zh**4))
    case(31)
        ! l=5, m=0, iang=31 (=l^2+l+1+m)
        ! Ylm= (Sqrt[11/Pi]*zh*(15 - 70*zh^2 + 63*zh^4))/16
        y=1.169503224534236E-1_r8*(zh*(15 - 70*zh**2 + 63*zh**4))
    case(32)
        ! l=5, m=1, iang=32 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[165/Pi]*xh*(1 - 14*zh^2 + 21*zh^4))/16
        y=-4.529466511956969E-1_r8*(xh*(1 - 14*zh**2 + 21*zh**4))
    case(33)
        ! l=5, m=2, iang=33 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*(xh^2 - yh^2)*zh*(-1 + 3*zh^2))/8
        y=2.396768392486662E0_r8*((xh**2 - yh**2)*zh*(-1 + 3*zh**2))
    case(34)
        ! l=5, m=3, iang=34 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 9*zh^2))/16
        y=-4.892382994352504E-1_r8*(xh*(xh**2 - 3*yh**2)*(-1 + 9*zh**2))
    case(35)
        ! l=5, m=4, iang=35 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh)/16
        y=2.075662314881041E0_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*zh)
    case(36)
        ! l=5, m=5, iang=36 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4))/16
        y=-6.563820568401701E-1_r8*(xh**5 - 10*xh**3*yh**2 + 5*xh*yh**4)
    case(37)
        ! l=6, m=-6, iang=37 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4))/16
        y=1.366368210383829E0_r8*(3*xh**5*yh - 10*xh**3*yh**3 + 3*xh*yh**5)
    case(38)
        ! l=6, m=-5, iang=38 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh)/16
        y=2.366619162231752E0_r8*(yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh)
    case(39)
        ! l=6, m=-4, iang=39 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*xh*yh*(xh^2 - yh^2)*(-1 + 11*zh^2))/8
        y=2.018259602914897E0_r8*(xh*yh*(xh**2 - yh**2)*(-1 + 11*zh**2))
    case(40)
        ! l=6, m=-3, iang=40 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(-3 + 11*zh^2))/16
        y=-9.212052595149235E-1_r8*(yh*(-3*xh**2 + yh**2)*zh*(-3 + 11*zh**2))
    case(41)
        ! l=6, m=-2, iang=41 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 33*zh^4))/16
        y=9.212052595149235E-1_r8*(xh*yh*(1 - 18*zh**2 + 33*zh**4))
    case(42)
        ! l=6, m=-1, iang=42 (=l^2+l+1+m)
        ! Ylm= (Sqrt[273/Pi]*yh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        y=5.826213625187314E-1_r8*(yh*zh*(5 - 30*zh**2 + 33*zh**4))
    case(43)
        ! l=6, m=0, iang=43 (=l^2+l+1+m)
        ! Ylm= (Sqrt[13/Pi]*(-5 + 105*zh^2 - 315*zh^4 + 231*zh^6))/32
        y=6.356920226762843E-2_r8*(-5 + 105*zh**2 - 315*zh**4 + 231*zh**6)
    case(44)
        ! l=6, m=1, iang=44 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[273/Pi]*xh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        y=-5.826213625187314E-1_r8*(xh*zh*(5 - 30*zh**2 + 33*zh**4))
    case(45)
        ! l=6, m=2, iang=45 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*(xh^2 - yh^2)*(1 - 18*zh^2 + 33*zh^4))/32
        y=4.606026297574617E-1_r8*((xh**2 - yh**2)*(1 - 18*zh**2 + 33*zh**4))
    case(46)
        ! l=6, m=3, iang=46 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(-3 + 11*zh^2))/16
        y=-9.212052595149235E-1_r8*(xh*(xh**2 - 3*yh**2)*zh*(-3 + 11*zh**2))
    case(47)
        ! l=6, m=4, iang=47 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 11*zh^2))/32
        y=5.045649007287242E-1_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 11*zh**2))
    case(48)
        ! l=6, m=5, iang=48 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh)/16
        y=-2.366619162231752E0_r8*(xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh)
    case(49)
        ! l=6, m=6, iang=49 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*(xh^6 - 15*xh^4*yh^2 + 15*xh^2*yh^4 - yh^6))/32
        y=6.831841051919143E-1_r8*(xh**6 - 15*xh**4*yh**2 + 15*xh**2*yh**4 - yh**6)
    end select
end subroutine

!> evaluate spherical harmonic functions + gradient
elemental subroutine spherical_harmonic_gradient(xh,yh,zh,ir,iang,y,gx,gy,gz)
    !> xhat, xh=x/r
    real(r8), intent(in) :: xh
    !> yhat, yh=y/r
    real(r8), intent(in) :: yh
    !> zhat, zh=z/r
    real(r8), intent(in) :: zh
    !> inverse norm of vector
    real(r8), intent(in) :: ir
    !> angular index, iang=l^2+l+m+1
    integer, intent(in) :: iang
    !> spherical harmonic
    real(r8), intent(out) :: y
    !> spherical harmonic gradient, x-direction
    real(r8), intent(out) :: gx
    !> spherical harmonic gradient, y-direction
    real(r8), intent(out) :: gy
    !> spherical harmonic gradient, z-direction
    real(r8), intent(out) :: gz

    select case(iang)
    case(1)
        ! l=0, m=0, iang=1 (=l^2+l+1+m)
        ! Ylm= 1/(2*Sqrt[Pi])
        ! dYlm/dx = 0
        ! dYlm/dy = 0
        ! dYlm/dz = 0
        y=2.820947917738781E-1_r8
        gx=0.0_r8
        gy=0.0_r8
        gz=0.0_r8
    case(2)
        ! l=1, m=-1, iang=2 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*yh)/2
        ! dYlm/dx = -(Sqrt[3/Pi]*xh*yh)/(2*r)
        ! dYlm/dy = -(Sqrt[3/Pi]*(-1 + yh^2))/(2*r)
        ! dYlm/dz = -(Sqrt[3/Pi]*yh*zh)/(2*r)
        y=4.886025119029199E-1_r8*(yh)
        gx=4.886025119029199E-1_r8*(-(ir*xh*yh))
        gy=4.886025119029199E-1_r8*(ir - ir*yh**2)
        gz=4.886025119029199E-1_r8*(-(ir*yh*zh))
    case(3)
        ! l=1, m=0, iang=3 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*zh)/2
        ! dYlm/dx = -(Sqrt[3/Pi]*xh*zh)/(2*r)
        ! dYlm/dy = -(Sqrt[3/Pi]*yh*zh)/(2*r)
        ! dYlm/dz = -(Sqrt[3/Pi]*(-1 + zh^2))/(2*r)
        y=4.886025119029199E-1_r8*(zh)
        gx=4.886025119029199E-1_r8*(-(ir*xh*zh))
        gy=4.886025119029199E-1_r8*(-(ir*yh*zh))
        gz=4.886025119029199E-1_r8*(ir - ir*zh**2)
    case(4)
        ! l=1, m=1, iang=4 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[3/Pi]*xh)/2
        ! dYlm/dx = (Sqrt[3/Pi]*(-1 + xh^2))/(2*r)
        ! dYlm/dy = (Sqrt[3/Pi]*xh*yh)/(2*r)
        ! dYlm/dz = (Sqrt[3/Pi]*xh*zh)/(2*r)
        y=-4.886025119029199E-1_r8*(xh)
        gx=-4.886025119029199E-1_r8*(ir - ir*xh**2)
        gy=-4.886025119029199E-1_r8*(-(ir*xh*yh))
        gz=-4.886025119029199E-1_r8*(-(ir*xh*zh))
    case(5)
        ! l=2, m=-2, iang=5 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*xh*yh)/2
        ! dYlm/dx = (Sqrt[15/Pi]*(1 - 2*xh^2)*yh)/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*xh*(1 - 2*yh^2))/(2*r)
        ! dYlm/dz = -((Sqrt[15/Pi]*xh*yh*zh)/r)
        y=1.092548430592079E0_r8*(xh*yh)
        gx=1.092548430592079E0_r8*(ir*(1 - 2*xh**2)*yh)
        gy=1.092548430592079E0_r8*(ir*xh*(1 - 2*yh**2))
        gz=1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
    case(6)
        ! l=2, m=-1, iang=6 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*yh*zh)/2
        ! dYlm/dx = -((Sqrt[15/Pi]*xh*yh*zh)/r)
        ! dYlm/dy = (Sqrt[15/Pi]*(1 - 2*yh^2)*zh)/(2*r)
        ! dYlm/dz = (Sqrt[15/Pi]*yh*(1 - 2*zh^2))/(2*r)
        y=1.092548430592079E0_r8*(yh*zh)
        gx=1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
        gy=1.092548430592079E0_r8*(ir*(1 - 2*yh**2)*zh)
        gz=1.092548430592079E0_r8*(ir*yh*(1 - 2*zh**2))
    case(7)
        ! l=2, m=0, iang=7 (=l^2+l+1+m)
        ! Ylm= (Sqrt[5/Pi]*(-1 + 3*zh^2))/4
        ! dYlm/dx = (Sqrt[5/Pi]*xh*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        ! dYlm/dy = (Sqrt[5/Pi]*yh*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        ! dYlm/dz = (Sqrt[5/Pi]*zh*(2 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        y=3.153915652525200E-1_r8*(-1 + 3*zh**2)
        gx=3.153915652525200E-1_r8*(2*ir*xh*(-1 + xh**2 + yh**2 - 2*zh**2))
        gy=3.153915652525200E-1_r8*(2*ir*yh*(-1 + xh**2 + yh**2 - 2*zh**2))
        gz=3.153915652525200E-1_r8*(2*ir*zh*(2 + xh**2 + yh**2 - 2*zh**2))
    case(8)
        ! l=2, m=1, iang=8 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[15/Pi]*xh*zh)/2
        ! dYlm/dx = (Sqrt[15/Pi]*(-1 + 2*xh^2)*zh)/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*xh*yh*zh)/r
        ! dYlm/dz = (Sqrt[15/Pi]*xh*(-1 + 2*zh^2))/(2*r)
        y=-1.092548430592079E0_r8*(xh*zh)
        gx=-1.092548430592079E0_r8*(ir*(1 - 2*xh**2)*zh)
        gy=-1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
        gz=-1.092548430592079E0_r8*(ir*xh*(1 - 2*zh**2))
    case(9)
        ! l=2, m=2, iang=9 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*(xh^2 - yh^2))/4
        ! dYlm/dx = (Sqrt[15/Pi]*xh*(1 - xh^2 + yh^2))/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*yh*(-1 - xh^2 + yh^2))/(2*r)
        ! dYlm/dz = (Sqrt[15/Pi]*(-xh^2 + yh^2)*zh)/(2*r)
        y=5.462742152960395E-1_r8*(xh**2 - yh**2)
        gx=5.462742152960395E-1_r8*(2*ir*xh*(1 - xh**2 + yh**2))
        gy=5.462742152960395E-1_r8*(2*ir*yh*(-1 - xh**2 + yh**2))
        gz=5.462742152960395E-1_r8*(2*ir*(-xh**2 + yh**2)*zh)
    case(10)
        ! l=3, m=-3, iang=10 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*xh*yh*(2 - 3*xh^2 + yh^2))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*(xh^2 - (1 + 3*xh^2)*yh^2 + yh^4))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh)/(4*r)
        y=-5.900435899266435E-1_r8*(-3*xh**2*yh + yh**3)
        gx=-5.900435899266435E-1_r8*(-3*ir*xh*yh*(2 - 3*xh**2 + yh**2))
        gy=-5.900435899266435E-1_r8*(-3*ir*(xh**2 - (1 + 3*xh**2)*yh**2 + yh**4))
        gz=-5.900435899266435E-1_r8*(-3*ir*yh*(-3*xh**2 + yh**2)*zh)
    case(11)
        ! l=3, m=-2, iang=11 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*xh*yh*zh)/2
        ! dYlm/dx = (Sqrt[105/Pi]*(1 - 3*xh^2)*yh*zh)/(2*r)
        ! dYlm/dy = (Sqrt[105/Pi]*xh*(1 - 3*yh^2)*zh)/(2*r)
        ! dYlm/dz = (Sqrt[105/Pi]*xh*yh*(1 - 3*zh^2))/(2*r)
        y=2.890611442640554E0_r8*(xh*yh*zh)
        gx=2.890611442640554E0_r8*(ir*(1 - 3*xh**2)*yh*zh)
        gy=2.890611442640554E0_r8*(ir*xh*(1 - 3*yh**2)*zh)
        gz=2.890611442640554E0_r8*(ir*xh*yh*(1 - 3*zh**2))
    case(12)
        ! l=3, m=-1, iang=12 (=l^2+l+1+m)
        ! Ylm= (Sqrt[21/(2*Pi)]*yh*(-1 + 5*zh^2))/4
        ! dYlm/dx = (Sqrt[21/(2*Pi)]*xh*yh*(-2 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        ! dYlm/dy = (Sqrt[21/(2*Pi)]*(xh^2*(-1 + 3*yh^2) + 4*zh^2 + 3*yh^2*(-1 + yh^2 - 4*zh^2)))/(4*r)
        ! dYlm/dz = (Sqrt[21/(2*Pi)]*yh*zh*(8 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        y=4.570457994644657E-1_r8*(yh*(-1 + 5*zh**2))
        gx=4.570457994644657E-1_r8*(ir*xh*yh*(-2 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gy=4.570457994644657E-1_r8*(ir*(xh**2*(-1 + 3*yh**2) + 4*zh**2 + 3*yh**2*(-1 + yh**2 - 4*zh**2)))
        gz=4.570457994644657E-1_r8*(ir*yh*zh*(8 + 3*xh**2 + 3*yh**2 - 12*zh**2))
    case(13)
        ! l=3, m=0, iang=13 (=l^2+l+1+m)
        ! Ylm= (Sqrt[7/Pi]*zh*(-3 + 5*zh^2))/4
        ! dYlm/dx = (3*Sqrt[7/Pi]*xh*zh*(1 - 5*zh^2))/(4*r)
        ! dYlm/dy = (3*Sqrt[7/Pi]*yh*zh*(1 - 5*zh^2))/(4*r)
        ! dYlm/dz = (-3*Sqrt[7/Pi]*(1 - 6*zh^2 + 5*zh^4))/(4*r)
        y=3.731763325901154E-1_r8*(zh*(-3 + 5*zh**2))
        gx=3.731763325901154E-1_r8*(3*ir*xh*zh*(1 - 5*zh**2))
        gy=3.731763325901154E-1_r8*(3*ir*yh*zh*(1 - 5*zh**2))
        gz=3.731763325901154E-1_r8*(-3*ir*(1 - 6*zh**2 + 5*zh**4))
    case(14)
        ! l=3, m=1, iang=14 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[21/(2*Pi)]*xh*(-1 + 5*zh^2))/4
        ! dYlm/dx = (Sqrt[21/(2*Pi)]*(yh^2 - 4*zh^2 - 3*xh^2*(-1 + xh^2 + yh^2 - 4*zh^2)))/(4*r)
        ! dYlm/dy = (Sqrt[21/(2*Pi)]*xh*yh*(2 - 3*xh^2 - 3*yh^2 + 12*zh^2))/(4*r)
        ! dYlm/dz = -(Sqrt[21/(2*Pi)]*xh*zh*(8 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        y=-4.570457994644657E-1_r8*(xh*(-1 + 5*zh**2))
        gx=-4.570457994644657E-1_r8*(-(ir*(yh**2 - 4*zh**2 - 3*xh**2*(-1 + xh**2 + yh**2 - 4*zh**2))))
        gy=-4.570457994644657E-1_r8*(ir*xh*yh*(-2 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gz=-4.570457994644657E-1_r8*(ir*xh*zh*(8 + 3*xh**2 + 3*yh**2 - 12*zh**2))
    case(15)
        ! l=3, m=2, iang=15 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*(xh^2 - yh^2)*zh)/4
        ! dYlm/dx = (Sqrt[105/Pi]*xh*(2 - 3*xh^2 + 3*yh^2)*zh)/(4*r)
        ! dYlm/dy = (Sqrt[105/Pi]*yh*(-2 - 3*xh^2 + 3*yh^2)*zh)/(4*r)
        ! dYlm/dz = -(Sqrt[105/Pi]*(xh - yh)*(xh + yh)*(-1 + 3*zh^2))/(4*r)
        y=1.445305721320277E0_r8*((xh**2 - yh**2)*zh)
        gx=1.445305721320277E0_r8*(ir*xh*(2 - 3*xh**2 + 3*yh**2)*zh)
        gy=1.445305721320277E0_r8*(ir*yh*(-2 - 3*xh**2 + 3*yh**2)*zh)
        gz=1.445305721320277E0_r8*(-(ir*(xh - yh)*(xh + yh)*(-1 + 3*zh**2)))
    case(16)
        ! l=3, m=3, iang=16 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*(xh^4 + yh^2 - xh^2*(1 + 3*yh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*xh*yh*(2 + xh^2 - 3*yh^2))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh)/(4*r)
        y=-5.900435899266435E-1_r8*(xh**3 - 3*xh*yh**2)
        gx=-5.900435899266435E-1_r8*(-3*ir*(xh**4 + yh**2 - xh**2*(1 + 3*yh**2)))
        gy=-5.900435899266435E-1_r8*(-3*ir*xh*yh*(2 + xh**2 - 3*yh**2))
        gz=-5.900435899266435E-1_r8*(-3*ir*xh*(xh**2 - 3*yh**2)*zh)
    case(17)
        ! l=4, m=-4, iang=17 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*xh*yh*(xh^2 - yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/Pi]*yh*(-4*xh^4 - yh^2 + xh^2*(3 + 4*yh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/Pi]*xh*(xh^2 - (3 + 4*xh^2)*yh^2 + 4*yh^4))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/Pi]*xh*yh*(-xh^2 + yh^2)*zh)/r
        y=2.503342941796705E0_r8*(xh*yh*(xh**2 - yh**2))
        gx=2.503342941796705E0_r8*(ir*yh*(-4*xh**4 - yh**2 + xh**2*(3 + 4*yh**2)))
        gy=2.503342941796705E0_r8*(ir*xh*(xh**2 - (3 + 4*xh**2)*yh**2 + 4*yh**4))
        gz=2.503342941796705E0_r8*(4*ir*xh*yh*(-xh**2 + yh**2)*zh)
    case(18)
        ! l=4, m=-3, iang=18 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 - 6*xh^2 + 2*yh^2)*zh)/(2*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*(3*xh^2 - 3*(1 + 4*xh^2)*yh^2 + 4*yh^4)*zh)/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 4*zh^2))/(4*r)
        y=-1.770130769779931E0_r8*(yh*(-3*xh**2 + yh**2)*zh)
        gx=-1.770130769779931E0_r8*(2*ir*xh*yh*(-3 + 6*xh**2 - 2*yh**2)*zh)
        gy=-1.770130769779931E0_r8*(-(ir*(3*xh**2 - 3*(1 + 4*xh**2)*yh**2 + 4*yh**4)*zh))
        gz=-1.770130769779931E0_r8*(-(ir*yh*(-3*xh**2 + yh**2)*(-1 + 4*zh**2)))
    case(19)
        ! l=4, m=-2, iang=19 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*xh*yh*(-1 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/Pi]*yh*(4*xh^4 - yh^2 + 6*zh^2 + xh^2*(-3 + 4*yh^2 - 24*zh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[5/Pi]*xh*(4*yh^4 + xh^2*(-1 + 4*yh^2) + 6*zh^2 - 3*yh^2*(1 + 8*zh^2)))/(4*r)
        ! dYlm/dz = (3*Sqrt[5/Pi]*xh*yh*zh*(3 + xh^2 + yh^2 - 6*zh^2))/r
        y=9.461746957575600E-1_r8*(xh*yh*(-1 + 7*zh**2))
        gx=9.461746957575600E-1_r8*(ir*yh*(4*xh**4 - yh**2 + 6*zh**2 + xh**2*(-3 + 4*yh**2 - 24*zh**2)))
        gy=9.461746957575600E-1_r8*(ir*xh*(4*yh**4 + xh**2*(-1 + 4*yh**2) + 6*zh**2 - 3*yh**2*(1 + 8*zh**2)))
        gz=9.461746957575600E-1_r8*(4*ir*xh*yh*zh*(3 + xh**2 + yh**2 - 6*zh**2))
    case(20)
        ! l=4, m=-1, iang=20 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/(2*Pi)]*yh*zh*(-3 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/(2*Pi)]*xh*yh*zh*(3 - 14*zh^2))/(2*r)
        ! dYlm/dy = (3*Sqrt[5/(2*Pi)]*zh*(-3 + 7*zh^2 + yh^2*(6 - 28*zh^2)))/(4*r)
        ! dYlm/dz = (-3*Sqrt[5/(2*Pi)]*yh*(3 - 27*zh^2 + 28*zh^4))/(4*r)
        y=6.690465435572892E-1_r8*(yh*zh*(-3 + 7*zh**2))
        gx=6.690465435572892E-1_r8*(2*ir*xh*yh*zh*(3 - 14*zh**2))
        gy=6.690465435572892E-1_r8*(ir*zh*(-3 + 7*zh**2 + yh**2*(6 - 28*zh**2)))
        gz=6.690465435572892E-1_r8*(-(ir*yh*(3 - 27*zh**2 + 28*zh**4)))
    case(21)
        ! l=4, m=0, iang=21 (=l^2+l+1+m)
        ! Ylm= (3*(3 - 30*zh^2 + 35*zh^4))/(16*Sqrt[Pi])
        ! dYlm/dx = (15*xh*zh^2*(3 - 7*zh^2))/(4*Sqrt[Pi]*r)
        ! dYlm/dy = (15*yh*zh^2*(3 - 7*zh^2))/(4*Sqrt[Pi]*r)
        ! dYlm/dz = (-15*zh*(3 - 10*zh^2 + 7*zh^4))/(4*Sqrt[Pi]*r)
        y=1.057855469152043E-1_r8*(3 - 30*zh**2 + 35*zh**4)
        gx=1.057855469152043E-1_r8*(20*ir*xh*zh**2*(3 - 7*zh**2))
        gy=1.057855469152043E-1_r8*(20*ir*yh*zh**2*(3 - 7*zh**2))
        gz=1.057855469152043E-1_r8*(-20*ir*zh*(3 - 10*zh**2 + 7*zh**4))
    case(22)
        ! l=4, m=1, iang=22 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[5/(2*Pi)]*xh*zh*(-3 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/(2*Pi)]*zh*(3 - 7*zh^2 + xh^2*(-6 + 28*zh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[5/(2*Pi)]*xh*yh*zh*(-3 + 14*zh^2))/(2*r)
        ! dYlm/dz = (3*Sqrt[5/(2*Pi)]*xh*(3 - 27*zh^2 + 28*zh^4))/(4*r)
        y=-6.690465435572892E-1_r8*(xh*zh*(-3 + 7*zh**2))
        gx=-6.690465435572892E-1_r8*(-(ir*zh*(3 - 7*zh**2 + xh**2*(-6 + 28*zh**2))))
        gy=-6.690465435572892E-1_r8*(-2*ir*xh*yh*zh*(-3 + 14*zh**2))
        gz=-6.690465435572892E-1_r8*(-(ir*xh*(3 - 27*zh**2 + 28*zh**4)))
    case(23)
        ! l=4, m=2, iang=23 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*(xh^2 - yh^2)*(-1 + 7*zh^2))/8
        ! dYlm/dx = (3*Sqrt[5/Pi]*xh*(-xh^2 + xh^4 - yh^4 + 3*(1 - 2*xh^2 + 2*yh^2)*zh^2))/(2*r)
        ! dYlm/dy = (3*Sqrt[5/Pi]*yh*(xh^4 + yh^2 - yh^4 - 3*(1 + 2*xh^2 - 2*yh^2)*zh^2))/(2*r)
        ! dYlm/dz = (3*Sqrt[5/Pi]*(xh - yh)*(xh + yh)*zh*(3 + xh^2 + yh^2 - 6*zh^2))/(2*r)
        y=4.730873478787800E-1_r8*((xh**2 - yh**2)*(-1 + 7*zh**2))
        gx=4.730873478787800E-1_r8*(4*ir*xh*(-xh**2 + xh**4 - yh**4 + 3*(1 - 2*xh**2 + 2*yh**2)*zh**2))
        gy=4.730873478787800E-1_r8*(4*ir*yh*(xh**4 + yh**2 - yh**4 - 3*(1 + 2*xh**2 - 2*yh**2)*zh**2))
        gz=4.730873478787800E-1_r8*(4*ir*(xh - yh)*(xh + yh)*zh*(3 + xh**2 + yh**2 - 6*zh**2))
    case(24)
        ! l=4, m=3, iang=24 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*(4*xh^4 + 3*yh^2 - 3*xh^2*(1 + 4*yh^2))*zh)/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 + 2*xh^2 - 6*yh^2)*zh)/(2*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 4*zh^2))/(4*r)
        y=-1.770130769779931E0_r8*(xh*(xh**2 - 3*yh**2)*zh)
        gx=-1.770130769779931E0_r8*(-(ir*(4*xh**4 + 3*yh**2 - 3*xh**2*(1 + 4*yh**2))*zh))
        gy=-1.770130769779931E0_r8*(-2*ir*xh*yh*(3 + 2*xh**2 - 6*yh**2)*zh)
        gz=-1.770130769779931E0_r8*(-(ir*xh*(xh**2 - 3*yh**2)*(-1 + 4*zh**2)))
    case(25)
        ! l=4, m=4, iang=25 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4))/16
        ! dYlm/dx = (-3*Sqrt[35/Pi]*xh*(xh^4 + 3*yh^2 + yh^4 - xh^2*(1 + 6*yh^2)))/(4*r)
        ! dYlm/dy = (-3*Sqrt[35/Pi]*yh*(xh^4 - yh^2 + yh^4 + xh^2*(3 - 6*yh^2)))/(4*r)
        ! dYlm/dz = (-3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh)/(4*r)
        y=6.258357354491761E-1_r8*(xh**4 - 6*xh**2*yh**2 + yh**4)
        gx=6.258357354491761E-1_r8*(-4*ir*xh*(xh**4 + 3*yh**2 + yh**4 - xh**2*(1 + 6*yh**2)))
        gy=6.258357354491761E-1_r8*(-4*ir*yh*(xh**4 - yh**2 + yh**4 + xh**2*(3 - 6*yh**2)))
        gz=6.258357354491761E-1_r8*(-4*ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*zh)
    case(26)
        ! l=5, m=-5, iang=26 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4))/16
        ! dYlm/dx = (-15*Sqrt[77/(2*Pi)]*xh*yh*(5*xh^4 + 4*yh^2 + yh^4 - 2*xh^2*(2 + 5*yh^2)))/(16*r)
        ! dYlm/dy = (-15*Sqrt[77/(2*Pi)]*(-yh^4 + yh^6 + 2*xh^2*yh^2*(3 - 5*yh^2) + xh^4*(-1 + 5*yh^2)))/(16*r)
        ! dYlm/dz = (-15*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh)/(16*r)
        y=6.563820568401701E-1_r8*(5*xh**4*yh - 10*xh**2*yh**3 + yh**5)
        gx=6.563820568401701E-1_r8*(-5*ir*xh*yh*(5*xh**4 + 4*yh**2 + yh**4 - 2*xh**2*(2 + 5*yh**2)))
        gy=6.563820568401701E-1_r8*(-5*ir*(-yh**4 + yh**6 + 2*xh**2*yh**2*(3 - 5*yh**2) + xh**4*(-1 + 5*yh**2)))
        gz=6.563820568401701E-1_r8*(-5*ir*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh)
    case(27)
        ! l=5, m=-4, iang=27 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*xh*yh*(xh^2 - yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[385/Pi]*yh*(-5*xh^4 - yh^2 + xh^2*(3 + 5*yh^2))*zh)/(4*r)
        ! dYlm/dy = (3*Sqrt[385/Pi]*xh*(xh^2 - (3 + 5*xh^2)*yh^2 + 5*yh^4)*zh)/(4*r)
        ! dYlm/dz = (-3*Sqrt[385/Pi]*xh*(xh - yh)*yh*(xh + yh)*(-1 + 5*zh^2))/(4*r)
        y=8.302649259524165E0_r8*(xh*yh*(xh**2 - yh**2)*zh)
        gx=8.302649259524165E0_r8*(ir*yh*(-5*xh**4 - yh**2 + xh**2*(3 + 5*yh**2))*zh)
        gy=8.302649259524165E0_r8*(ir*xh*(xh**2 - (3 + 5*xh**2)*yh**2 + 5*yh**4)*zh)
        gz=8.302649259524165E0_r8*(-(ir*xh*(xh - yh)*yh*(xh + yh)*(-1 + 5*zh**2)))
    case(28)
        ! l=5, m=-3, iang=28 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 9*zh^2))/16
        ! dYlm/dx = (Sqrt[385/(2*Pi)]*xh*yh*(15*xh^4 - 5*yh^4 + 48*zh^2 + 2*xh^2*(-6 + 5*yh^2 - 60*zh^2) + yh^2*(-4 + 40*zh^2)))/(16*r)
        ! dYlm/dy = (Sqrt[385/(2*Pi)]*(-5*yh^6 + 3*xh^4*(-1 + 5*yh^2) - 24*yh^2*zh^2 + 5*yh^4*(1 + 8*zh^2) + 2*xh^2*(5*yh^4 + 12*zh^2 - 3*yh^2*(1 + 20*zh^2))))/(16*r)
        ! dYlm/dz = -(Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(16 + 5*xh^2 + 5*yh^2 - 40*zh^2))/(16*r)
        y=-4.892382994352504E-1_r8*(yh*(-3*xh**2 + yh**2)*(-1 + 9*zh**2))
        gx=-4.892382994352504E-1_r8*(-(ir*xh*yh*(15*xh**4 - 5*yh**4 + 48*zh**2 + 2*xh**2*(-6 + 5*yh**2 - 60*zh**2) + yh**2*(-4 + 40*zh**2))))
        gy=-4.892382994352504E-1_r8*(-(ir*(-5*yh**6 + 3*xh**4*(-1 + 5*yh**2) - 24*yh**2*zh**2 + 5*yh**4*(1 + 8*zh**2) + 2*xh**2*(5*yh**4 + 12*zh**2 - 3*yh**2*(1 + 20*zh**2)))))
        gz=-4.892382994352504E-1_r8*(ir*yh*(-3*xh**2 + yh**2)*zh*(16 + 5*xh**2 + 5*yh**2 - 40*zh**2))
    case(29)
        ! l=5, m=-2, iang=29 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*xh*yh*zh*(-1 + 3*zh^2))/4
        ! dYlm/dx = (Sqrt[1155/Pi]*yh*zh*(5*xh^4 - yh^2 + 2*zh^2 + xh^2*(-3 + 5*yh^2 - 10*zh^2)))/(4*r)
        ! dYlm/dy = (Sqrt[1155/Pi]*xh*zh*(5*yh^4 + xh^2*(-1 + 5*yh^2) + 2*zh^2 - yh^2*(3 + 10*zh^2)))/(4*r)
        ! dYlm/dz = (Sqrt[1155/Pi]*xh*yh*(-xh^2 - yh^2 + (6 + 5*xh^2 + 5*yh^2)*zh^2 - 10*zh^4))/(4*r)
        y=4.793536784973324E0_r8*(xh*yh*zh*(-1 + 3*zh**2))
        gx=4.793536784973324E0_r8*(ir*yh*zh*(5*xh**4 - yh**2 + 2*zh**2 + xh**2*(-3 + 5*yh**2 - 10*zh**2)))
        gy=4.793536784973324E0_r8*(ir*xh*zh*(5*yh**4 + xh**2*(-1 + 5*yh**2) + 2*zh**2 - yh**2*(3 + 10*zh**2)))
        gz=4.793536784973324E0_r8*(ir*xh*yh*(-xh**2 - yh**2 + (6 + 5*xh**2 + 5*yh**2)*zh**2 - 10*zh**4))
    case(30)
        ! l=5, m=-1, iang=30 (=l^2+l+1+m)
        ! Ylm= (Sqrt[165/Pi]*yh*(1 - 14*zh^2 + 21*zh^4))/16
        ! dYlm/dx = -(Sqrt[165/Pi]*xh*yh*(1 - 42*zh^2 + 105*zh^4))/(16*r)
        ! dYlm/dy = -(Sqrt[165/Pi]*(-1 + yh^2 + 14*(1 - 3*yh^2)*zh^2 + 21*(-1 + 5*yh^2)*zh^4))/(16*r)
        ! dYlm/dz = (Sqrt[165/Pi]*yh*zh*(-29 + 21*zh^2*(6 - 5*zh^2)))/(16*r)
        y=4.529466511956969E-1_r8*(yh*(1 - 14*zh**2 + 21*zh**4))
        gx=4.529466511956969E-1_r8*(-(ir*xh*yh*(1 - 42*zh**2 + 105*zh**4)))
        gy=4.529466511956969E-1_r8*(-(ir*(-1 + yh**2 + 14*(1 - 3*yh**2)*zh**2 + 21*(-1 + 5*yh**2)*zh**4)))
        gz=4.529466511956969E-1_r8*(ir*yh*zh*(-29 + 21*zh**2*(6 - 5*zh**2)))
    case(31)
        ! l=5, m=0, iang=31 (=l^2+l+1+m)
        ! Ylm= (Sqrt[11/Pi]*zh*(15 - 70*zh^2 + 63*zh^4))/16
        ! dYlm/dx = (-15*Sqrt[11/Pi]*xh*zh*(1 - 14*zh^2 + 21*zh^4))/(16*r)
        ! dYlm/dy = (-15*Sqrt[11/Pi]*yh*zh*(1 - 14*zh^2 + 21*zh^4))/(16*r)
        ! dYlm/dz = (15*Sqrt[11/Pi]*(1 - 15*zh^2 + 35*zh^4 - 21*zh^6))/(16*r)
        y=1.169503224534236E-1_r8*(zh*(15 - 70*zh**2 + 63*zh**4))
        gx=1.169503224534236E-1_r8*(-15*ir*xh*zh*(1 - 14*zh**2 + 21*zh**4))
        gy=1.169503224534236E-1_r8*(-15*ir*yh*zh*(1 - 14*zh**2 + 21*zh**4))
        gz=1.169503224534236E-1_r8*(15*ir*(1 - 15*zh**2 + 35*zh**4 - 21*zh**6))
    case(32)
        ! l=5, m=1, iang=32 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[165/Pi]*xh*(1 - 14*zh^2 + 21*zh^4))/16
        ! dYlm/dx = (Sqrt[165/Pi]*(-1 + xh^2 + 14*(1 - 3*xh^2)*zh^2 + 21*(-1 + 5*xh^2)*zh^4))/(16*r)
        ! dYlm/dy = (Sqrt[165/Pi]*xh*yh*(1 - 42*zh^2 + 105*zh^4))/(16*r)
        ! dYlm/dz = (Sqrt[165/Pi]*xh*zh*(29 + 21*zh^2*(-6 + 5*zh^2)))/(16*r)
        y=-4.529466511956969E-1_r8*(xh*(1 - 14*zh**2 + 21*zh**4))
        gx=-4.529466511956969E-1_r8*(-(ir*(-1 + xh**2 + 14*(1 - 3*xh**2)*zh**2 + 21*(-1 + 5*xh**2)*zh**4)))
        gy=-4.529466511956969E-1_r8*(-(ir*xh*yh*(1 - 42*zh**2 + 105*zh**4)))
        gz=-4.529466511956969E-1_r8*(-(ir*xh*zh*(29 + 21*zh**2*(-6 + 5*zh**2))))
    case(33)
        ! l=5, m=2, iang=33 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*(xh^2 - yh^2)*zh*(-1 + 3*zh^2))/8
        ! dYlm/dx = (Sqrt[1155/Pi]*xh*zh*(-4*xh^2 + 5*xh^4 - 5*yh^4 + 2*(2 - 5*xh^2 + 5*yh^2)*zh^2))/(8*r)
        ! dYlm/dy = (Sqrt[1155/Pi]*yh*zh*(5*xh^4 + 4*yh^2 - 5*yh^4 - 2*(2 + 5*xh^2 - 5*yh^2)*zh^2))/(8*r)
        ! dYlm/dz = (Sqrt[1155/Pi]*(xh - yh)*(xh + yh)*(-xh^2 - yh^2 + (6 + 5*xh^2 + 5*yh^2)*zh^2 - 10*zh^4))/(8*r)
        y=2.396768392486662E0_r8*((xh**2 - yh**2)*zh*(-1 + 3*zh**2))
        gx=2.396768392486662E0_r8*(ir*xh*zh*(-4*xh**2 + 5*xh**4 - 5*yh**4 + 2*(2 - 5*xh**2 + 5*yh**2)*zh**2))
        gy=2.396768392486662E0_r8*(ir*yh*zh*(5*xh**4 + 4*yh**2 - 5*yh**4 - 2*(2 + 5*xh**2 - 5*yh**2)*zh**2))
        gz=2.396768392486662E0_r8*(ir*(xh - yh)*(xh + yh)*(-xh**2 - yh**2 + (6 + 5*xh**2 + 5*yh**2)*zh**2 - 10*zh**4))
    case(34)
        ! l=5, m=3, iang=34 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 9*zh^2))/16
        ! dYlm/dx = -(Sqrt[385/(2*Pi)]*(5*xh^6 - 5*xh^4*(1 + 2*yh^2 + 8*zh^2) + 3*(yh^4 - 8*yh^2*zh^2) + 3*xh^2*(-5*yh^4 + 8*zh^2 + yh^2*(2 + 40*zh^2))))/(16*r)
        ! dYlm/dy = (Sqrt[385/(2*Pi)]*xh*yh*(-5*xh^4 + 48*zh^2 + 3*yh^2*(-4 + 5*yh^2 - 40*zh^2) + 2*xh^2*(-2 + 5*yh^2 + 20*zh^2)))/(16*r)
        ! dYlm/dz = -(Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(16 + 5*xh^2 + 5*yh^2 - 40*zh^2))/(16*r)
        y=-4.892382994352504E-1_r8*(xh*(xh**2 - 3*yh**2)*(-1 + 9*zh**2))
        gx=-4.892382994352504E-1_r8*(ir*(5*xh**6 - 5*xh**4*(1 + 2*yh**2 + 8*zh**2) + 3*(yh**4 - 8*yh**2*zh**2) + 3*xh**2*(-5*yh**4 + 8*zh**2 + yh**2*(2 + 40*zh**2))))
        gy=-4.892382994352504E-1_r8*(-(ir*xh*yh*(-5*xh**4 + 48*zh**2 + 3*yh**2*(-4 + 5*yh**2 - 40*zh**2) + 2*xh**2*(-2 + 5*yh**2 + 20*zh**2))))
        gz=-4.892382994352504E-1_r8*(ir*xh*(xh**2 - 3*yh**2)*zh*(16 + 5*xh**2 + 5*yh**2 - 40*zh**2))
    case(35)
        ! l=5, m=4, iang=35 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh)/16
        ! dYlm/dx = (-3*Sqrt[385/Pi]*xh*(xh^2*(-4 + 5*xh^2) + 6*(2 - 5*xh^2)*yh^2 + 5*yh^4)*zh)/(16*r)
        ! dYlm/dy = (-3*Sqrt[385/Pi]*yh*(5*xh^4 - 4*yh^2 + 5*yh^4 + 6*xh^2*(2 - 5*yh^2))*zh)/(16*r)
        ! dYlm/dz = (-3*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 5*zh^2))/(16*r)
        y=2.075662314881041E0_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*zh)
        gx=2.075662314881041E0_r8*(-(ir*xh*(xh**2*(-4 + 5*xh**2) + 6*(2 - 5*xh**2)*yh**2 + 5*yh**4)*zh))
        gy=2.075662314881041E0_r8*(-(ir*yh*(5*xh**4 - 4*yh**2 + 5*yh**4 + 6*xh**2*(2 - 5*yh**2))*zh))
        gz=2.075662314881041E0_r8*(-(ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 5*zh**2)))
    case(36)
        ! l=5, m=5, iang=36 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4))/16
        ! dYlm/dx = (15*Sqrt[77/(2*Pi)]*(xh^6 - yh^4 + xh^2*yh^2*(6 + 5*yh^2) - xh^4*(1 + 10*yh^2)))/(16*r)
        ! dYlm/dy = (15*Sqrt[77/(2*Pi)]*xh*yh*(xh^2*(4 + xh^2) - 2*(2 + 5*xh^2)*yh^2 + 5*yh^4))/(16*r)
        ! dYlm/dz = (15*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh)/(16*r)
        y=-6.563820568401701E-1_r8*(xh**5 - 10*xh**3*yh**2 + 5*xh*yh**4)
        gx=-6.563820568401701E-1_r8*(-5*ir*(xh**6 - yh**4 + xh**2*yh**2*(6 + 5*yh**2) - xh**4*(1 + 10*yh**2)))
        gy=-6.563820568401701E-1_r8*(-5*ir*xh*yh*(xh**2*(4 + xh**2) - 2*(2 + 5*xh**2)*yh**2 + 5*yh**4))
        gz=-6.563820568401701E-1_r8*(-5*ir*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh)
    case(37)
        ! l=6, m=-6, iang=37 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4))/16
        ! dYlm/dx = (3*Sqrt[3003/(2*Pi)]*yh*(-6*xh^6 + yh^4 - 2*xh^2*yh^2*(5 + 3*yh^2) + 5*xh^4*(1 + 4*yh^2)))/(16*r)
        ! dYlm/dy = (3*Sqrt[3003/(2*Pi)]*xh*(xh^4 - 2*xh^2*(5 + 3*xh^2)*yh^2 + 5*(1 + 4*xh^2)*yh^4 - 6*yh^6))/(16*r)
        ! dYlm/dz = (-3*Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4)*zh)/(8*r)
        y=1.366368210383829E0_r8*(3*xh**5*yh - 10*xh**3*yh**3 + 3*xh*yh**5)
        gx=1.366368210383829E0_r8*(3*ir*yh*(-6*xh**6 + yh**4 - 2*xh**2*yh**2*(5 + 3*yh**2) + 5*xh**4*(1 + 4*yh**2)))
        gy=1.366368210383829E0_r8*(3*ir*xh*(xh**4 - 2*xh**2*(5 + 3*xh**2)*yh**2 + 5*(1 + 4*xh**2)*yh**4 - 6*yh**6))
        gz=1.366368210383829E0_r8*(-6*ir*xh*yh*(3*xh**4 - 10*xh**2*yh**2 + 3*yh**4)*zh)
    case(38)
        ! l=6, m=-5, iang=38 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh)/16
        ! dYlm/dx = (-3*Sqrt[1001/(2*Pi)]*xh*yh*(15*xh^4 + 10*yh^2 + 3*yh^4 - 10*xh^2*(1 + 3*yh^2))*zh)/(8*r)
        ! dYlm/dy = (-3*Sqrt[1001/(2*Pi)]*(-5*xh^4 + 30*xh^2*(1 + xh^2)*yh^2 - 5*(1 + 12*xh^2)*yh^4 + 6*yh^6)*zh)/(16*r)
        ! dYlm/dz = (-3*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*(-1 + 6*zh^2))/(16*r)
        y=2.366619162231752E0_r8*(yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh)
        gx=2.366619162231752E0_r8*(-2*ir*xh*yh*(15*xh**4 + 10*yh**2 + 3*yh**4 - 10*xh**2*(1 + 3*yh**2))*zh)
        gy=2.366619162231752E0_r8*(-(ir*(-5*xh**4 + 30*xh**2*(1 + xh**2)*yh**2 - 5*(1 + 12*xh**2)*yh**4 + 6*yh**6)*zh))
        gz=2.366619162231752E0_r8*(-(ir*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*(-1 + 6*zh**2)))
    case(39)
        ! l=6, m=-4, iang=39 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*xh*yh*(xh^2 - yh^2)*(-1 + 11*zh^2))/8
        ! dYlm/dx = (3*Sqrt[91/Pi]*yh*(-5*xh^4 + 6*xh^6 + yh^4 - 6*xh^2*yh^4 - 10*(6*xh^4 + yh^2 - 3*xh^2*(1 + 2*yh^2))*zh^2))/(8*r)
        ! dYlm/dy = (3*Sqrt[91/Pi]*xh*(-xh^4 + 6*xh^4*yh^2 + 5*yh^4 - 6*yh^6 + 10*(xh^2 - 3*(1 + 2*xh^2)*yh^2 + 6*yh^4)*zh^2))/(8*r)
        ! dYlm/dz = (3*Sqrt[91/Pi]*xh*(xh - yh)*yh*(xh + yh)*zh*(10 + 3*xh^2 + 3*yh^2 - 30*zh^2))/(4*r)
        y=2.018259602914897E0_r8*(xh*yh*(xh**2 - yh**2)*(-1 + 11*zh**2))
        gx=2.018259602914897E0_r8*(ir*yh*(-5*xh**4 + 6*xh**6 + yh**4 - 6*xh**2*yh**4 - 10*(6*xh**4 + yh**2 - 3*xh**2*(1 + 2*yh**2))*zh**2))
        gy=2.018259602914897E0_r8*(ir*xh*(-xh**4 + 6*xh**4*yh**2 + 5*yh**4 - 6*yh**6 + 10*(xh**2 - 3*(1 + 2*xh**2)*yh**2 + 6*yh**4)*zh**2))
        gz=2.018259602914897E0_r8*(2*ir*xh*(xh - yh)*yh*(xh + yh)*zh*(10 + 3*xh**2 + 3*yh**2 - 30*zh**2))
    case(40)
        ! l=6, m=-3, iang=40 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(-3 + 11*zh^2))/16
        ! dYlm/dx = (-3*Sqrt[1365/(2*Pi)]*xh*yh*zh*(3 - 6*xh^2 + 2*yh^2 - 11*(1 - 3*xh^2 + yh^2)*zh^2))/(8*r)
        ! dYlm/dy = (3*Sqrt[1365/(2*Pi)]*zh*(xh^2*(-3 + 11*zh^2 + yh^2*(12 - 66*zh^2)) + yh^2*(3 - 11*zh^2 + yh^2*(-4 + 22*zh^2))))/(16*r)
        ! dYlm/dz = (3*Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(1 - 15*zh^2 + 22*zh^4))/(16*r)
        y=-9.212052595149235E-1_r8*(yh*(-3*xh**2 + yh**2)*zh*(-3 + 11*zh**2))
        gx=-9.212052595149235E-1_r8*(6*ir*xh*yh*zh*(3 - 6*xh**2 + 2*yh**2 - 11*(1 - 3*xh**2 + yh**2)*zh**2))
        gy=-9.212052595149235E-1_r8*(-3*ir*zh*(xh**2*(-3 + 11*zh**2 + yh**2*(12 - 66*zh**2)) + yh**2*(3 - 11*zh**2 + yh**2*(-4 + 22*zh**2))))
        gz=-9.212052595149235E-1_r8*(-3*ir*yh*(-3*xh**2 + yh**2)*(1 - 15*zh**2 + 22*zh**4))
    case(41)
        ! l=6, m=-2, iang=41 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 33*zh^4))/16
        ! dYlm/dx = (Sqrt[1365/(2*Pi)]*yh*(1 - 18*zh^2 + 33*zh^4 - 2*xh^2*(1 - 36*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dy = (Sqrt[1365/(2*Pi)]*xh*(1 - 18*zh^2 + 33*zh^4 - 2*yh^2*(1 - 36*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dz = -(Sqrt[1365/(2*Pi)]*xh*yh*zh*(19 - 102*zh^2 + 99*zh^4))/(8*r)
        y=9.212052595149235E-1_r8*(xh*yh*(1 - 18*zh**2 + 33*zh**4))
        gx=9.212052595149235E-1_r8*(ir*yh*(1 - 18*zh**2 + 33*zh**4 - 2*xh**2*(1 - 36*zh**2 + 99*zh**4)))
        gy=9.212052595149235E-1_r8*(ir*xh*(1 - 18*zh**2 + 33*zh**4 - 2*yh**2*(1 - 36*zh**2 + 99*zh**4)))
        gz=9.212052595149235E-1_r8*(-2*ir*xh*yh*zh*(19 - 102*zh**2 + 99*zh**4))
    case(42)
        ! l=6, m=-1, iang=42 (=l^2+l+1+m)
        ! Ylm= (Sqrt[273/Pi]*yh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        ! dYlm/dx = -(Sqrt[273/Pi]*xh*yh*zh*(5 - 60*zh^2 + 99*zh^4))/(8*r)
        ! dYlm/dy = (Sqrt[273/Pi]*zh*(5 - 30*zh^2 + 33*zh^4 - 2*yh^2*(5 - 60*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dz = (Sqrt[273/Pi]*yh*(5 - 100*zh^2 + 285*zh^4 - 198*zh^6))/(16*r)
        y=5.826213625187314E-1_r8*(yh*zh*(5 - 30*zh**2 + 33*zh**4))
        gx=5.826213625187314E-1_r8*(-2*ir*xh*yh*zh*(5 - 60*zh**2 + 99*zh**4))
        gy=5.826213625187314E-1_r8*(ir*zh*(5 - 30*zh**2 + 33*zh**4 - 2*yh**2*(5 - 60*zh**2 + 99*zh**4)))
        gz=5.826213625187314E-1_r8*(ir*yh*(5 - 100*zh**2 + 285*zh**4 - 198*zh**6))
    case(43)
        ! l=6, m=0, iang=43 (=l^2+l+1+m)
        ! Ylm= (Sqrt[13/Pi]*(-5 + 105*zh^2 - 315*zh^4 + 231*zh^6))/32
        ! dYlm/dx = (-21*Sqrt[13/Pi]*xh*zh^2*(5 - 30*zh^2 + 33*zh^4))/(16*r)
        ! dYlm/dy = (-21*Sqrt[13/Pi]*yh*zh^2*(5 - 30*zh^2 + 33*zh^4))/(16*r)
        ! dYlm/dz = (-21*Sqrt[13/Pi]*zh*(-5 + 35*zh^2 - 63*zh^4 + 33*zh^6))/(16*r)
        y=6.356920226762843E-2_r8*(-5 + 105*zh**2 - 315*zh**4 + 231*zh**6)
        gx=6.356920226762843E-2_r8*(-42*ir*xh*zh**2*(5 - 30*zh**2 + 33*zh**4))
        gy=6.356920226762843E-2_r8*(-42*ir*yh*zh**2*(5 - 30*zh**2 + 33*zh**4))
        gz=6.356920226762843E-2_r8*(-42*ir*zh*(-5 + 35*zh**2 - 63*zh**4 + 33*zh**6))
    case(44)
        ! l=6, m=1, iang=44 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[273/Pi]*xh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        ! dYlm/dx = (Sqrt[273/Pi]*zh*(-5 + 30*zh^2 - 33*zh^4 + 2*xh^2*(5 - 60*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dy = (Sqrt[273/Pi]*xh*yh*zh*(5 - 60*zh^2 + 99*zh^4))/(8*r)
        ! dYlm/dz = (Sqrt[273/Pi]*xh*(-5 + 100*zh^2 - 285*zh^4 + 198*zh^6))/(16*r)
        y=-5.826213625187314E-1_r8*(xh*zh*(5 - 30*zh**2 + 33*zh**4))
        gx=-5.826213625187314E-1_r8*(-(ir*zh*(-5 + 30*zh**2 - 33*zh**4 + 2*xh**2*(5 - 60*zh**2 + 99*zh**4))))
        gy=-5.826213625187314E-1_r8*(-2*ir*xh*yh*zh*(5 - 60*zh**2 + 99*zh**4))
        gz=-5.826213625187314E-1_r8*(-(ir*xh*(-5 + 100*zh**2 - 285*zh**4 + 198*zh**6)))
    case(45)
        ! l=6, m=2, iang=45 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*(xh^2 - yh^2)*(1 - 18*zh^2 + 33*zh^4))/32
        ! dYlm/dx = (Sqrt[1365/(2*Pi)]*xh*(1 - xh^2 + yh^2 + 18*(-1 + 2*xh^2 - 2*yh^2)*zh^2 + 33*(1 - 3*xh^2 + 3*yh^2)*zh^4))/(16*r)
        ! dYlm/dy = (Sqrt[1365/(2*Pi)]*yh*(-1 - xh^2 + yh^2 + 18*(1 + 2*xh^2 - 2*yh^2)*zh^2 - 33*(1 + 3*xh^2 - 3*yh^2)*zh^4))/(16*r)
        ! dYlm/dz = -(Sqrt[1365/(2*Pi)]*(xh - yh)*(xh + yh)*zh*(19 - 102*zh^2 + 99*zh^4))/(16*r)
        y=4.606026297574617E-1_r8*((xh**2 - yh**2)*(1 - 18*zh**2 + 33*zh**4))
        gx=4.606026297574617E-1_r8*(2*ir*xh*(1 - xh**2 + yh**2 + 18*(-1 + 2*xh**2 - 2*yh**2)*zh**2 + 33*(1 - 3*xh**2 + 3*yh**2)*zh**4))
        gy=4.606026297574617E-1_r8*(2*ir*yh*(-1 - xh**2 + yh**2 + 18*(1 + 2*xh**2 - 2*yh**2)*zh**2 - 33*(1 + 3*xh**2 - 3*yh**2)*zh**4))
        gz=4.606026297574617E-1_r8*(-2*ir*(xh - yh)*(xh + yh)*zh*(19 - 102*zh**2 + 99*zh**4))
    case(46)
        ! l=6, m=3, iang=46 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(-3 + 11*zh^2))/16
        ! dYlm/dx = (3*Sqrt[1365/(2*Pi)]*zh*(yh^2*(-3 + 11*zh^2) + xh^4*(-4 + 22*zh^2) + xh^2*(3 + 12*yh^2 - 11*(1 + 6*yh^2)*zh^2)))/(16*r)
        ! dYlm/dy = (3*Sqrt[1365/(2*Pi)]*xh*yh*zh*(-3 - 2*xh^2 + 6*yh^2 + 11*(1 + xh^2 - 3*yh^2)*zh^2))/(8*r)
        ! dYlm/dz = (3*Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(1 - 15*zh^2 + 22*zh^4))/(16*r)
        y=-9.212052595149235E-1_r8*(xh*(xh**2 - 3*yh**2)*zh*(-3 + 11*zh**2))
        gx=-9.212052595149235E-1_r8*(-3*ir*zh*(yh**2*(-3 + 11*zh**2) + xh**4*(-4 + 22*zh**2) + xh**2*(3 + 12*yh**2 - 11*(1 + 6*yh**2)*zh**2)))
        gy=-9.212052595149235E-1_r8*(-6*ir*xh*yh*zh*(-3 - 2*xh**2 + 6*yh**2 + 11*(1 + xh**2 - 3*yh**2)*zh**2))
        gz=-9.212052595149235E-1_r8*(-3*ir*xh*(xh**2 - 3*yh**2)*(1 - 15*zh**2 + 22*zh**4))
    case(47)
        ! l=6, m=4, iang=47 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 11*zh^2))/32
        ! dYlm/dx = (3*Sqrt[91/Pi]*xh*(3*xh^6 + yh^4*(5 + 3*yh^2) - 30*yh^2*(2 + yh^2)*zh^2 - 3*xh^4*(1 + 5*yh^2 + 10*zh^2) + 5*xh^2*(-3*yh^4 + 4*zh^2 + yh^2*(2 + 36*zh^2))))/(16*r)
        ! dYlm/dy = (3*Sqrt[91/Pi]*yh*(3*xh^6 + 3*yh^6 + 20*yh^2*zh^2 - 5*xh^4*(-1 + 3*yh^2 + 6*zh^2) - 3*yh^4*(1 + 10*zh^2) + 5*xh^2*(-3*yh^4 - 12*zh^2 + yh^2*(2 + 36*zh^2))))/(16*r)
        ! dYlm/dz = (3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh*(10 + 3*xh^2 + 3*yh^2 - 30*zh^2))/(16*r)
        y=5.045649007287242E-1_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 11*zh**2))
        gx=5.045649007287242E-1_r8*(2*ir*xh*(3*xh**6 + yh**4*(5 + 3*yh**2) - 30*yh**2*(2 + yh**2)*zh**2 - 3*xh**4*(1 + 5*yh**2 + 10*zh**2) + 5*xh**2*(-3*yh**4 + 4*zh**2 + yh**2*(2 + 36*zh**2))))
        gy=5.045649007287242E-1_r8*(2*ir*yh*(3*xh**6 + 3*yh**6 + 20*yh**2*zh**2 - 5*xh**4*(-1 + 3*yh**2 + 6*zh**2) - 3*yh**4*(1 + 10*zh**2) + 5*xh**2*(-3*yh**4 - 12*zh**2 + yh**2*(2 + 36*zh**2))))
        gz=5.045649007287242E-1_r8*(2*ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*zh*(10 + 3*xh**2 + 3*yh**2 - 30*zh**2))
    case(48)
        ! l=6, m=5, iang=48 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh)/16
        ! dYlm/dx = (3*Sqrt[1001/(2*Pi)]*(6*xh^6 - 5*yh^4 - 5*xh^4*(1 + 12*yh^2) + 30*xh^2*(yh^2 + yh^4))*zh)/(16*r)
        ! dYlm/dy = (3*Sqrt[1001/(2*Pi)]*xh*yh*(3*xh^4 + xh^2*(10 - 30*yh^2) + 5*yh^2*(-2 + 3*yh^2))*zh)/(8*r)
        ! dYlm/dz = (3*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*(-1 + 6*zh^2))/(16*r)
        y=-2.366619162231752E0_r8*(xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh)
        gx=-2.366619162231752E0_r8*(-(ir*(6*xh**6 - 5*yh**4 - 5*xh**4*(1 + 12*yh**2) + 30*xh**2*(yh**2 + yh**4))*zh))
        gy=-2.366619162231752E0_r8*(-2*ir*xh*yh*(3*xh**4 + xh**2*(10 - 30*yh**2) + 5*yh**2*(-2 + 3*yh**2))*zh)
        gz=-2.366619162231752E0_r8*(-(ir*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*(-1 + 6*zh**2)))
    case(49)
        ! l=6, m=6, iang=49 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*(xh^6 - 15*xh^4*yh^2 + 15*xh^2*yh^4 - yh^6))/32
        ! dYlm/dx = (-3*Sqrt[3003/(2*Pi)]*(xh^7 - xh*yh^4*(5 + yh^2) + 5*xh^3*yh^2*(2 + 3*yh^2) - xh^5*(1 + 15*yh^2)))/(16*r)
        ! dYlm/dy = (-3*Sqrt[3003/(2*Pi)]*yh*(xh^6 + yh^4 - yh^6 + xh^4*(5 - 15*yh^2) + 5*xh^2*yh^2*(-2 + 3*yh^2)))/(16*r)
        ! dYlm/dz = (3*Sqrt[3003/(2*Pi)]*(-xh^6 + 15*xh^4*yh^2 - 15*xh^2*yh^4 + yh^6)*zh)/(16*r)
        y=6.831841051919143E-1_r8*(xh**6 - 15*xh**4*yh**2 + 15*xh**2*yh**4 - yh**6)
        gx=6.831841051919143E-1_r8*(-6*ir*(xh**7 - xh*yh**4*(5 + yh**2) + 5*xh**3*yh**2*(2 + 3*yh**2) - xh**5*(1 + 15*yh**2)))
        gy=6.831841051919143E-1_r8*(-6*ir*yh*(xh**6 + yh**4 - yh**6 + xh**4*(5 - 15*yh**2) + 5*xh**2*yh**2*(-2 + 3*yh**2)))
        gz=6.831841051919143E-1_r8*(6*ir*(-xh**6 + 15*xh**4*yh**2 - 15*xh**2*yh**4 + yh**6)*zh)
    end select
end subroutine

!> evaluate spherical harmonic functions + gradient
elemental subroutine spherical_harmonic_gradient_hessian(xh,yh,zh,ir,iang,y,gx,gy,gz,gxx,gxy,gxz,gyy,gyz,gzz)
    !> xhat, xh=x/r
    real(r8), intent(in) :: xh
    !> yhat, yh=y/r
    real(r8), intent(in) :: yh
    !> zhat, zh=z/r
    real(r8), intent(in) :: zh
    !> inverse norm of vector
    real(r8), intent(in) :: ir
    !> angular index, iang=l^2+l+m+1
    integer, intent(in) :: iang
    !> spherical harmonic
    real(r8), intent(out) :: y
    !> spherical harmonic gradient, x-direction
    real(r8), intent(out) :: gx
    !> spherical harmonic gradient, y-direction
    real(r8), intent(out) :: gy
    !> spherical harmonic gradient, z-direction
    real(r8), intent(out) :: gz
    !> spherical harmonic hessian, xx-direction
    real(r8), intent(out) :: gxx
    !> spherical harmonic hessian, xy-direction
    real(r8), intent(out) :: gxy
    !> spherical harmonic hessian, xz-direction
    real(r8), intent(out) :: gxz
    !> spherical harmonic hessian, yy-direction
    real(r8), intent(out) :: gyy
    !> spherical harmonic hessian, yz-direction
    real(r8), intent(out) :: gyz
    !> spherical harmonic hessian, zz-direction
    real(r8), intent(out) :: gzz

    select case(iang)
    case(1)
        ! l=0, m=0, iang=1 (=l^2+l+1+m)
        ! Ylm= 1/(2*Sqrt[Pi])
        ! dYlm/dx = 0
        ! dYlm/dy = 0
        ! dYlm/dz = 0
        ! d^2Ylm/dxdx = 0
        ! d^2Ylm/dxdy = 0
        ! d^2Ylm/dxdz = 0
        ! d^2Ylm/dydx = 0
        ! d^2Ylm/dydy = 0
        ! d^2Ylm/dydz = 0
        ! d^2Ylm/dzdx = 0
        ! d^2Ylm/dzdy = 0
        ! d^2Ylm/dzdz = 0
        y=2.820947917738781E-1_r8
        gx=0.0_r8
        gy=0.0_r8
        gz=0.0_r8
        gxx=0.0_r8
        gxy=0.0_r8
        gxz=0.0_r8
        gyy=0.0_r8
        gyz=0.0_r8
        gzz=0.0_r8
    case(2)
        ! l=1, m=-1, iang=2 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*yh)/2
        ! dYlm/dx = -(Sqrt[3/Pi]*xh*yh)/(2*r)
        ! dYlm/dy = -(Sqrt[3/Pi]*(-1 + yh^2))/(2*r)
        ! dYlm/dz = -(Sqrt[3/Pi]*yh*zh)/(2*r)
        ! d^2Ylm/dxdx = (Sqrt[3/Pi]*(-1 + 3*xh^2)*yh)/(2*r^2)
        ! d^2Ylm/dxdy = (Sqrt[3/Pi]*xh*(-1 + 3*yh^2))/(2*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dydx = (Sqrt[3/Pi]*xh*(-1 + 3*yh^2))/(2*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[3/Pi]*yh*(-1 + yh^2))/(2*r^2)
        ! d^2Ylm/dydz = (Sqrt[3/Pi]*(-1 + 3*yh^2)*zh)/(2*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dzdy = (Sqrt[3/Pi]*(-1 + 3*yh^2)*zh)/(2*r^2)
        ! d^2Ylm/dzdz = (Sqrt[3/Pi]*yh*(-1 + 3*zh^2))/(2*r^2)
        y=4.886025119029199E-1_r8*(yh)
        gx=4.886025119029199E-1_r8*(-(ir*xh*yh))
        gy=4.886025119029199E-1_r8*(ir - ir*yh**2)
        gz=4.886025119029199E-1_r8*(-(ir*yh*zh))
        gxx=4.886025119029199E-1_r8*(ir**2*(-1 + 3*xh**2)*yh)
        gxy=4.886025119029199E-1_r8*(ir**2*xh*(-1 + 3*yh**2))
        gxz=4.886025119029199E-1_r8*(3*ir**2*xh*yh*zh)
        gyy=4.886025119029199E-1_r8*(3*ir**2*yh*(-1 + yh**2))
        gyz=4.886025119029199E-1_r8*(ir**2*(-1 + 3*yh**2)*zh)
        gzz=4.886025119029199E-1_r8*(ir**2*yh*(-1 + 3*zh**2))
    case(3)
        ! l=1, m=0, iang=3 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3/Pi]*zh)/2
        ! dYlm/dx = -(Sqrt[3/Pi]*xh*zh)/(2*r)
        ! dYlm/dy = -(Sqrt[3/Pi]*yh*zh)/(2*r)
        ! dYlm/dz = -(Sqrt[3/Pi]*(-1 + zh^2))/(2*r)
        ! d^2Ylm/dxdx = (Sqrt[3/Pi]*(-1 + 3*xh^2)*zh)/(2*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dxdz = (Sqrt[3/Pi]*xh*(-1 + 3*zh^2))/(2*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dydy = (Sqrt[3/Pi]*(-1 + 3*yh^2)*zh)/(2*r^2)
        ! d^2Ylm/dydz = (Sqrt[3/Pi]*yh*(-1 + 3*zh^2))/(2*r^2)
        ! d^2Ylm/dzdx = (Sqrt[3/Pi]*xh*(-1 + 3*zh^2))/(2*r^2)
        ! d^2Ylm/dzdy = (Sqrt[3/Pi]*yh*(-1 + 3*zh^2))/(2*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[3/Pi]*zh*(-1 + zh^2))/(2*r^2)
        y=4.886025119029199E-1_r8*(zh)
        gx=4.886025119029199E-1_r8*(-(ir*xh*zh))
        gy=4.886025119029199E-1_r8*(-(ir*yh*zh))
        gz=4.886025119029199E-1_r8*(ir - ir*zh**2)
        gxx=4.886025119029199E-1_r8*(ir**2*(-1 + 3*xh**2)*zh)
        gxy=4.886025119029199E-1_r8*(3*ir**2*xh*yh*zh)
        gxz=4.886025119029199E-1_r8*(ir**2*xh*(-1 + 3*zh**2))
        gyy=4.886025119029199E-1_r8*(ir**2*(-1 + 3*yh**2)*zh)
        gyz=4.886025119029199E-1_r8*(ir**2*yh*(-1 + 3*zh**2))
        gzz=4.886025119029199E-1_r8*(3*ir**2*zh*(-1 + zh**2))
    case(4)
        ! l=1, m=1, iang=4 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[3/Pi]*xh)/2
        ! dYlm/dx = (Sqrt[3/Pi]*(-1 + xh^2))/(2*r)
        ! dYlm/dy = (Sqrt[3/Pi]*xh*yh)/(2*r)
        ! dYlm/dz = (Sqrt[3/Pi]*xh*zh)/(2*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[3/Pi]*xh*(-1 + xh^2))/(2*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[3/Pi]*(-1 + 3*xh^2)*yh)/(2*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[3/Pi]*(-1 + 3*xh^2)*zh)/(2*r^2)
        ! d^2Ylm/dydx = -(Sqrt[3/Pi]*(-1 + 3*xh^2)*yh)/(2*r^2)
        ! d^2Ylm/dydy = -(Sqrt[3/Pi]*xh*(-1 + 3*yh^2))/(2*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[3/Pi]*(-1 + 3*xh^2)*zh)/(2*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[3/Pi]*xh*yh*zh)/(2*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[3/Pi]*xh*(-1 + 3*zh^2))/(2*r^2)
        y=-4.886025119029199E-1_r8*(xh)
        gx=-4.886025119029199E-1_r8*(ir - ir*xh**2)
        gy=-4.886025119029199E-1_r8*(-(ir*xh*yh))
        gz=-4.886025119029199E-1_r8*(-(ir*xh*zh))
        gxx=-4.886025119029199E-1_r8*(3*ir**2*xh*(-1 + xh**2))
        gxy=-4.886025119029199E-1_r8*(ir**2*(-1 + 3*xh**2)*yh)
        gxz=-4.886025119029199E-1_r8*(ir**2*(-1 + 3*xh**2)*zh)
        gyy=-4.886025119029199E-1_r8*(ir**2*xh*(-1 + 3*yh**2))
        gyz=-4.886025119029199E-1_r8*(3*ir**2*xh*yh*zh)
        gzz=-4.886025119029199E-1_r8*(ir**2*xh*(-1 + 3*zh**2))
    case(5)
        ! l=2, m=-2, iang=5 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*xh*yh)/2
        ! dYlm/dx = (Sqrt[15/Pi]*(1 - 2*xh^2)*yh)/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*xh*(1 - 2*yh^2))/(2*r)
        ! dYlm/dz = -((Sqrt[15/Pi]*xh*yh*zh)/r)
        ! d^2Ylm/dxdx = (Sqrt[15/Pi]*xh*(-3 + 4*xh^2)*yh)/r^2
        ! d^2Ylm/dxdy = (Sqrt[15/Pi]*(1 - 2*yh^2 + xh^2*(-2 + 8*yh^2)))/(2*r^2)
        ! d^2Ylm/dxdz = (Sqrt[15/Pi]*(-1 + 4*xh^2)*yh*zh)/r^2
        ! d^2Ylm/dydx = (Sqrt[15/Pi]*(1 - 2*yh^2 + xh^2*(-2 + 8*yh^2)))/(2*r^2)
        ! d^2Ylm/dydy = (Sqrt[15/Pi]*xh*yh*(-3 + 4*yh^2))/r^2
        ! d^2Ylm/dydz = (Sqrt[15/Pi]*xh*(-1 + 4*yh^2)*zh)/r^2
        ! d^2Ylm/dzdx = (Sqrt[15/Pi]*(-1 + 4*xh^2)*yh*zh)/r^2
        ! d^2Ylm/dzdy = (Sqrt[15/Pi]*xh*(-1 + 4*yh^2)*zh)/r^2
        ! d^2Ylm/dzdz = (Sqrt[15/Pi]*xh*yh*(-1 + 4*zh^2))/r^2
        y=1.092548430592079E0_r8*(xh*yh)
        gx=1.092548430592079E0_r8*(ir*(1 - 2*xh**2)*yh)
        gy=1.092548430592079E0_r8*(ir*xh*(1 - 2*yh**2))
        gz=1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
        gxx=1.092548430592079E0_r8*(2*ir**2*xh*(-3 + 4*xh**2)*yh)
        gxy=1.092548430592079E0_r8*(ir**2*(1 - 2*yh**2 + xh**2*(-2 + 8*yh**2)))
        gxz=1.092548430592079E0_r8*(2*ir**2*(-1 + 4*xh**2)*yh*zh)
        gyy=1.092548430592079E0_r8*(2*ir**2*xh*yh*(-3 + 4*yh**2))
        gyz=1.092548430592079E0_r8*(2*ir**2*xh*(-1 + 4*yh**2)*zh)
        gzz=1.092548430592079E0_r8*(2*ir**2*xh*yh*(-1 + 4*zh**2))
    case(6)
        ! l=2, m=-1, iang=6 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*yh*zh)/2
        ! dYlm/dx = -((Sqrt[15/Pi]*xh*yh*zh)/r)
        ! dYlm/dy = (Sqrt[15/Pi]*(1 - 2*yh^2)*zh)/(2*r)
        ! dYlm/dz = (Sqrt[15/Pi]*yh*(1 - 2*zh^2))/(2*r)
        ! d^2Ylm/dxdx = (Sqrt[15/Pi]*(-1 + 4*xh^2)*yh*zh)/r^2
        ! d^2Ylm/dxdy = (Sqrt[15/Pi]*xh*(-1 + 4*yh^2)*zh)/r^2
        ! d^2Ylm/dxdz = (Sqrt[15/Pi]*xh*yh*(-1 + 4*zh^2))/r^2
        ! d^2Ylm/dydx = (Sqrt[15/Pi]*xh*(-1 + 4*yh^2)*zh)/r^2
        ! d^2Ylm/dydy = (Sqrt[15/Pi]*yh*(-3 + 4*yh^2)*zh)/r^2
        ! d^2Ylm/dydz = (Sqrt[15/Pi]*(1 - 2*zh^2 + yh^2*(-2 + 8*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdx = (Sqrt[15/Pi]*xh*yh*(-1 + 4*zh^2))/r^2
        ! d^2Ylm/dzdy = (Sqrt[15/Pi]*(1 - 2*zh^2 + yh^2*(-2 + 8*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdz = (Sqrt[15/Pi]*yh*zh*(-3 + 4*zh^2))/r^2
        y=1.092548430592079E0_r8*(yh*zh)
        gx=1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
        gy=1.092548430592079E0_r8*(ir*(1 - 2*yh**2)*zh)
        gz=1.092548430592079E0_r8*(ir*yh*(1 - 2*zh**2))
        gxx=1.092548430592079E0_r8*(2*ir**2*(-1 + 4*xh**2)*yh*zh)
        gxy=1.092548430592079E0_r8*(2*ir**2*xh*(-1 + 4*yh**2)*zh)
        gxz=1.092548430592079E0_r8*(2*ir**2*xh*yh*(-1 + 4*zh**2))
        gyy=1.092548430592079E0_r8*(2*ir**2*yh*(-3 + 4*yh**2)*zh)
        gyz=1.092548430592079E0_r8*(ir**2*(1 - 2*zh**2 + yh**2*(-2 + 8*zh**2)))
        gzz=1.092548430592079E0_r8*(2*ir**2*yh*zh*(-3 + 4*zh**2))
    case(7)
        ! l=2, m=0, iang=7 (=l^2+l+1+m)
        ! Ylm= (Sqrt[5/Pi]*(-1 + 3*zh^2))/4
        ! dYlm/dx = (Sqrt[5/Pi]*xh*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        ! dYlm/dy = (Sqrt[5/Pi]*yh*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        ! dYlm/dz = (Sqrt[5/Pi]*zh*(2 + xh^2 + yh^2 - 2*zh^2))/(2*r)
        ! d^2Ylm/dxdx = -(Sqrt[5/Pi]*(-1 + 4*xh^2)*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r^2)
        ! d^2Ylm/dxdy = (-2*Sqrt[5/Pi]*xh*yh*(-1 + xh^2 + yh^2 - 2*zh^2))/r^2
        ! d^2Ylm/dxdz = -((Sqrt[5/Pi]*xh*zh*(1 + 2*xh^2 + 2*yh^2 - 4*zh^2))/r^2)
        ! d^2Ylm/dydx = (-2*Sqrt[5/Pi]*xh*yh*(-1 + xh^2 + yh^2 - 2*zh^2))/r^2
        ! d^2Ylm/dydy = -(Sqrt[5/Pi]*(-1 + 4*yh^2)*(-1 + xh^2 + yh^2 - 2*zh^2))/(2*r^2)
        ! d^2Ylm/dydz = -((Sqrt[5/Pi]*yh*zh*(1 + 2*xh^2 + 2*yh^2 - 4*zh^2))/r^2)
        ! d^2Ylm/dzdx = -((Sqrt[5/Pi]*xh*zh*(1 + 2*xh^2 + 2*yh^2 - 4*zh^2))/r^2)
        ! d^2Ylm/dzdy = -((Sqrt[5/Pi]*yh*zh*(1 + 2*xh^2 + 2*yh^2 - 4*zh^2))/r^2)
        ! d^2Ylm/dzdz = -(Sqrt[5/Pi]*(2 + xh^2 + yh^2 - 2*zh^2)*(-1 + 4*zh^2))/(2*r^2)
        y=3.153915652525200E-1_r8*(-1 + 3*zh**2)
        gx=3.153915652525200E-1_r8*(2*ir*xh*(-1 + xh**2 + yh**2 - 2*zh**2))
        gy=3.153915652525200E-1_r8*(2*ir*yh*(-1 + xh**2 + yh**2 - 2*zh**2))
        gz=3.153915652525200E-1_r8*(2*ir*zh*(2 + xh**2 + yh**2 - 2*zh**2))
        gxx=3.153915652525200E-1_r8*(-2*ir**2*(-1 + 4*xh**2)*(-1 + xh**2 + yh**2 - 2*zh**2))
        gxy=3.153915652525200E-1_r8*(-8*ir**2*xh*yh*(-1 + xh**2 + yh**2 - 2*zh**2))
        gxz=3.153915652525200E-1_r8*(-4*ir**2*xh*zh*(1 + 2*xh**2 + 2*yh**2 - 4*zh**2))
        gyy=3.153915652525200E-1_r8*(-2*ir**2*(-1 + 4*yh**2)*(-1 + xh**2 + yh**2 - 2*zh**2))
        gyz=3.153915652525200E-1_r8*(-4*ir**2*yh*zh*(1 + 2*xh**2 + 2*yh**2 - 4*zh**2))
        gzz=3.153915652525200E-1_r8*(-2*ir**2*(2 + xh**2 + yh**2 - 2*zh**2)*(-1 + 4*zh**2))
    case(8)
        ! l=2, m=1, iang=8 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[15/Pi]*xh*zh)/2
        ! dYlm/dx = (Sqrt[15/Pi]*(-1 + 2*xh^2)*zh)/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*xh*yh*zh)/r
        ! dYlm/dz = (Sqrt[15/Pi]*xh*(-1 + 2*zh^2))/(2*r)
        ! d^2Ylm/dxdx = (Sqrt[15/Pi]*xh*(3 - 4*xh^2)*zh)/r^2
        ! d^2Ylm/dxdy = (Sqrt[15/Pi]*(1 - 4*xh^2)*yh*zh)/r^2
        ! d^2Ylm/dxdz = (Sqrt[15/Pi]*(-1 + 2*zh^2 + xh^2*(2 - 8*zh^2)))/(2*r^2)
        ! d^2Ylm/dydx = (Sqrt[15/Pi]*(1 - 4*xh^2)*yh*zh)/r^2
        ! d^2Ylm/dydy = (Sqrt[15/Pi]*xh*(1 - 4*yh^2)*zh)/r^2
        ! d^2Ylm/dydz = (Sqrt[15/Pi]*xh*yh*(1 - 4*zh^2))/r^2
        ! d^2Ylm/dzdx = (Sqrt[15/Pi]*(-1 + 2*zh^2 + xh^2*(2 - 8*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdy = (Sqrt[15/Pi]*xh*yh*(1 - 4*zh^2))/r^2
        ! d^2Ylm/dzdz = (Sqrt[15/Pi]*xh*zh*(3 - 4*zh^2))/r^2
        y=-1.092548430592079E0_r8*(xh*zh)
        gx=-1.092548430592079E0_r8*(ir*(1 - 2*xh**2)*zh)
        gy=-1.092548430592079E0_r8*(-2*ir*xh*yh*zh)
        gz=-1.092548430592079E0_r8*(ir*xh*(1 - 2*zh**2))
        gxx=-1.092548430592079E0_r8*(2*ir**2*xh*(-3 + 4*xh**2)*zh)
        gxy=-1.092548430592079E0_r8*(2*ir**2*(-1 + 4*xh**2)*yh*zh)
        gxz=-1.092548430592079E0_r8*(ir**2*(1 - 2*zh**2 + xh**2*(-2 + 8*zh**2)))
        gyy=-1.092548430592079E0_r8*(2*ir**2*xh*(-1 + 4*yh**2)*zh)
        gyz=-1.092548430592079E0_r8*(2*ir**2*xh*yh*(-1 + 4*zh**2))
        gzz=-1.092548430592079E0_r8*(2*ir**2*xh*zh*(-3 + 4*zh**2))
    case(9)
        ! l=2, m=2, iang=9 (=l^2+l+1+m)
        ! Ylm= (Sqrt[15/Pi]*(xh^2 - yh^2))/4
        ! dYlm/dx = (Sqrt[15/Pi]*xh*(1 - xh^2 + yh^2))/(2*r)
        ! dYlm/dy = (Sqrt[15/Pi]*yh*(-1 - xh^2 + yh^2))/(2*r)
        ! dYlm/dz = (Sqrt[15/Pi]*(-xh^2 + yh^2)*zh)/(2*r)
        ! d^2Ylm/dxdx = (Sqrt[15/Pi]*(-1 + 4*xh^2)*(-1 + xh^2 - yh^2))/(2*r^2)
        ! d^2Ylm/dxdy = (2*Sqrt[15/Pi]*xh*yh*(xh^2 - yh^2))/r^2
        ! d^2Ylm/dxdz = (Sqrt[15/Pi]*xh*(-1 + 2*xh^2 - 2*yh^2)*zh)/r^2
        ! d^2Ylm/dydx = (2*Sqrt[15/Pi]*xh*yh*(xh^2 - yh^2))/r^2
        ! d^2Ylm/dydy = (Sqrt[15/Pi]*(1 + xh^2 - yh^2)*(-1 + 4*yh^2))/(2*r^2)
        ! d^2Ylm/dydz = -((Sqrt[15/Pi]*yh*(-1 - 2*xh^2 + 2*yh^2)*zh)/r^2)
        ! d^2Ylm/dzdx = (Sqrt[15/Pi]*xh*(-1 + 2*xh^2 - 2*yh^2)*zh)/r^2
        ! d^2Ylm/dzdy = -((Sqrt[15/Pi]*yh*(-1 - 2*xh^2 + 2*yh^2)*zh)/r^2)
        ! d^2Ylm/dzdz = (Sqrt[15/Pi]*(xh^2 - yh^2)*(-1 + 4*zh^2))/(2*r^2)
        y=5.462742152960395E-1_r8*(xh**2 - yh**2)
        gx=5.462742152960395E-1_r8*(2*ir*xh*(1 - xh**2 + yh**2))
        gy=5.462742152960395E-1_r8*(2*ir*yh*(-1 - xh**2 + yh**2))
        gz=5.462742152960395E-1_r8*(2*ir*(-xh**2 + yh**2)*zh)
        gxx=5.462742152960395E-1_r8*(2*ir**2*(-1 + 4*xh**2)*(-1 + xh**2 - yh**2))
        gxy=5.462742152960395E-1_r8*(8*ir**2*xh*yh*(xh**2 - yh**2))
        gxz=5.462742152960395E-1_r8*(4*ir**2*xh*(-1 + 2*xh**2 - 2*yh**2)*zh)
        gyy=5.462742152960395E-1_r8*(2*ir**2*(1 + xh**2 - yh**2)*(-1 + 4*yh**2))
        gyz=5.462742152960395E-1_r8*(-4*ir**2*yh*(-1 - 2*xh**2 + 2*yh**2)*zh)
        gzz=5.462742152960395E-1_r8*(2*ir**2*(xh**2 - yh**2)*(-1 + 4*zh**2))
    case(10)
        ! l=3, m=-3, iang=10 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*xh*yh*(2 - 3*xh^2 + yh^2))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*(xh^2 - (1 + 3*xh^2)*yh^2 + yh^4))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh)/(4*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[35/(2*Pi)]*yh*(-2 - 15*xh^4 - yh^2 + 5*xh^2*(3 + yh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[35/(2*Pi)]*xh*(2 - 3*yh^2 - 5*yh^4 + 3*xh^2*(-1 + 5*yh^2)))/(4*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[35/(2*Pi)]*xh*yh*(-6 + 15*xh^2 - 5*yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[35/(2*Pi)]*xh*(2 - 3*yh^2 - 5*yh^4 + 3*xh^2*(-1 + 5*yh^2)))/(4*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[35/(2*Pi)]*yh*(2 - 7*yh^2 + 5*yh^4 + xh^2*(9 - 15*yh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[35/(2*Pi)]*(3*yh^2 - 5*yh^4 + 3*xh^2*(-1 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[35/(2*Pi)]*xh*yh*(-6 + 15*xh^2 - 5*yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[35/(2*Pi)]*(3*yh^2 - 5*yh^4 + 3*xh^2*(-1 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 5*zh^2))/(4*r^2)
        y=-5.900435899266435E-1_r8*(-3*xh**2*yh + yh**3)
        gx=-5.900435899266435E-1_r8*(-3*ir*xh*yh*(2 - 3*xh**2 + yh**2))
        gy=-5.900435899266435E-1_r8*(-3*ir*(xh**2 - (1 + 3*xh**2)*yh**2 + yh**4))
        gz=-5.900435899266435E-1_r8*(-3*ir*yh*(-3*xh**2 + yh**2)*zh)
        gxx=-5.900435899266435E-1_r8*(3*ir**2*yh*(-2 - 15*xh**4 - yh**2 + 5*xh**2*(3 + yh**2)))
        gxy=-5.900435899266435E-1_r8*(-3*ir**2*xh*(2 - 3*yh**2 - 5*yh**4 + 3*xh**2*(-1 + 5*yh**2)))
        gxz=-5.900435899266435E-1_r8*(-3*ir**2*xh*yh*(-6 + 15*xh**2 - 5*yh**2)*zh)
        gyy=-5.900435899266435E-1_r8*(3*ir**2*yh*(2 - 7*yh**2 + 5*yh**4 + xh**2*(9 - 15*yh**2)))
        gyz=-5.900435899266435E-1_r8*(-3*ir**2*(3*yh**2 - 5*yh**4 + 3*xh**2*(-1 + 5*yh**2))*zh)
        gzz=-5.900435899266435E-1_r8*(3*ir**2*yh*(-3*xh**2 + yh**2)*(-1 + 5*zh**2))
    case(11)
        ! l=3, m=-2, iang=11 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*xh*yh*zh)/2
        ! dYlm/dx = (Sqrt[105/Pi]*(1 - 3*xh^2)*yh*zh)/(2*r)
        ! dYlm/dy = (Sqrt[105/Pi]*xh*(1 - 3*yh^2)*zh)/(2*r)
        ! dYlm/dz = (Sqrt[105/Pi]*xh*yh*(1 - 3*zh^2))/(2*r)
        ! d^2Ylm/dxdx = (3*Sqrt[105/Pi]*xh*(-3 + 5*xh^2)*yh*zh)/(2*r^2)
        ! d^2Ylm/dxdy = (Sqrt[105/Pi]*(1 - 3*yh^2 + 3*xh^2*(-1 + 5*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dxdz = (Sqrt[105/Pi]*yh*(1 - 3*zh^2 + 3*xh^2*(-1 + 5*zh^2)))/(2*r^2)
        ! d^2Ylm/dydx = (Sqrt[105/Pi]*(1 - 3*yh^2 + 3*xh^2*(-1 + 5*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[105/Pi]*xh*yh*(-3 + 5*yh^2)*zh)/(2*r^2)
        ! d^2Ylm/dydz = (Sqrt[105/Pi]*xh*(1 - 3*zh^2 + 3*yh^2*(-1 + 5*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdx = (Sqrt[105/Pi]*yh*(1 - 3*zh^2 + 3*xh^2*(-1 + 5*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdy = (Sqrt[105/Pi]*xh*(1 - 3*zh^2 + 3*yh^2*(-1 + 5*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[105/Pi]*xh*yh*zh*(-3 + 5*zh^2))/(2*r^2)
        y=2.890611442640554E0_r8*(xh*yh*zh)
        gx=2.890611442640554E0_r8*(ir*(1 - 3*xh**2)*yh*zh)
        gy=2.890611442640554E0_r8*(ir*xh*(1 - 3*yh**2)*zh)
        gz=2.890611442640554E0_r8*(ir*xh*yh*(1 - 3*zh**2))
        gxx=2.890611442640554E0_r8*(3*ir**2*xh*(-3 + 5*xh**2)*yh*zh)
        gxy=2.890611442640554E0_r8*(ir**2*(1 - 3*yh**2 + 3*xh**2*(-1 + 5*yh**2))*zh)
        gxz=2.890611442640554E0_r8*(ir**2*yh*(1 - 3*zh**2 + 3*xh**2*(-1 + 5*zh**2)))
        gyy=2.890611442640554E0_r8*(3*ir**2*xh*yh*(-3 + 5*yh**2)*zh)
        gyz=2.890611442640554E0_r8*(ir**2*xh*(1 - 3*zh**2 + 3*yh**2*(-1 + 5*zh**2)))
        gzz=2.890611442640554E0_r8*(3*ir**2*xh*yh*zh*(-3 + 5*zh**2))
    case(12)
        ! l=3, m=-1, iang=12 (=l^2+l+1+m)
        ! Ylm= (Sqrt[21/(2*Pi)]*yh*(-1 + 5*zh^2))/4
        ! dYlm/dx = (Sqrt[21/(2*Pi)]*xh*yh*(-2 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        ! dYlm/dy = (Sqrt[21/(2*Pi)]*(xh^2*(-1 + 3*yh^2) + 4*zh^2 + 3*yh^2*(-1 + yh^2 - 4*zh^2)))/(4*r)
        ! dYlm/dz = (Sqrt[21/(2*Pi)]*yh*zh*(8 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        ! d^2Ylm/dxdx = -(Sqrt[21/(2*Pi)]*yh*(2 + 15*xh^4 - 3*yh^2 + 12*zh^2 + 15*xh^2*(-1 + yh^2 - 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[21/(2*Pi)]*xh*(2 + 15*yh^4 + 3*xh^2*(-1 + 5*yh^2) + 12*zh^2 - 15*yh^2*(1 + 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[21/(2*Pi)]*xh*yh*zh*(6 + 5*xh^2 + 5*yh^2 - 20*zh^2))/(4*r^2)
        ! d^2Ylm/dydx = -(Sqrt[21/(2*Pi)]*xh*(2 + 15*yh^4 + 3*xh^2*(-1 + 5*yh^2) + 12*zh^2 - 15*yh^2*(1 + 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[21/(2*Pi)]*yh*(2 + 5*yh^4 + xh^2*(-3 + 5*yh^2) + 12*zh^2 - yh^2*(7 + 20*zh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (Sqrt[21/(2*Pi)]*zh*(8 - 15*yh^4 + xh^2*(3 - 15*yh^2) - 12*zh^2 + 15*yh^2*(-1 + 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdx = (-3*Sqrt[21/(2*Pi)]*xh*yh*zh*(6 + 5*xh^2 + 5*yh^2 - 20*zh^2))/(4*r^2)
        ! d^2Ylm/dzdy = (Sqrt[21/(2*Pi)]*zh*(8 - 15*yh^4 + xh^2*(3 - 15*yh^2) - 12*zh^2 + 15*yh^2*(-1 + 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[21/(2*Pi)]*yh*(-8 + 60*zh^2 - 60*zh^4 + 3*xh^2*(-1 + 5*zh^2) + 3*yh^2*(-1 + 5*zh^2)))/(4*r^2)
        y=4.570457994644657E-1_r8*(yh*(-1 + 5*zh**2))
        gx=4.570457994644657E-1_r8*(ir*xh*yh*(-2 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gy=4.570457994644657E-1_r8*(ir*(xh**2*(-1 + 3*yh**2) + 4*zh**2 + 3*yh**2*(-1 + yh**2 - 4*zh**2)))
        gz=4.570457994644657E-1_r8*(ir*yh*zh*(8 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gxx=4.570457994644657E-1_r8*(-(ir**2*yh*(2 + 15*xh**4 - 3*yh**2 + 12*zh**2 + 15*xh**2*(-1 + yh**2 - 4*zh**2))))
        gxy=4.570457994644657E-1_r8*(-(ir**2*xh*(2 + 15*yh**4 + 3*xh**2*(-1 + 5*yh**2) + 12*zh**2 - 15*yh**2*(1 + 4*zh**2))))
        gxz=4.570457994644657E-1_r8*(-3*ir**2*xh*yh*zh*(6 + 5*xh**2 + 5*yh**2 - 20*zh**2))
        gyy=4.570457994644657E-1_r8*(-3*ir**2*yh*(2 + 5*yh**4 + xh**2*(-3 + 5*yh**2) + 12*zh**2 - yh**2*(7 + 20*zh**2)))
        gyz=4.570457994644657E-1_r8*(ir**2*zh*(8 - 15*yh**4 + xh**2*(3 - 15*yh**2) - 12*zh**2 + 15*yh**2*(-1 + 4*zh**2)))
        gzz=4.570457994644657E-1_r8*(-(ir**2*yh*(-8 + 60*zh**2 - 60*zh**4 + 3*xh**2*(-1 + 5*zh**2) + 3*yh**2*(-1 + 5*zh**2))))
    case(13)
        ! l=3, m=0, iang=13 (=l^2+l+1+m)
        ! Ylm= (Sqrt[7/Pi]*zh*(-3 + 5*zh^2))/4
        ! dYlm/dx = (3*Sqrt[7/Pi]*xh*zh*(1 - 5*zh^2))/(4*r)
        ! dYlm/dy = (3*Sqrt[7/Pi]*yh*zh*(1 - 5*zh^2))/(4*r)
        ! dYlm/dz = (-3*Sqrt[7/Pi]*(1 - 6*zh^2 + 5*zh^4))/(4*r)
        ! d^2Ylm/dxdx = (3*Sqrt[7/Pi]*zh*(1 - 5*zh^2 + xh^2*(-3 + 25*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[7/Pi]*xh*yh*zh*(-3 + 25*zh^2))/(4*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[7/Pi]*xh*(1 - 18*zh^2 + 25*zh^4))/(4*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[7/Pi]*xh*yh*zh*(-3 + 25*zh^2))/(4*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[7/Pi]*zh*(1 - 5*zh^2 + yh^2*(-3 + 25*zh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[7/Pi]*yh*(1 - 18*zh^2 + 25*zh^4))/(4*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[7/Pi]*xh*(1 - 18*zh^2 + 25*zh^4))/(4*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[7/Pi]*yh*(1 - 18*zh^2 + 25*zh^4))/(4*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[7/Pi]*zh*(13 - 38*zh^2 + 25*zh^4))/(4*r^2)
        y=3.731763325901154E-1_r8*(zh*(-3 + 5*zh**2))
        gx=3.731763325901154E-1_r8*(3*ir*xh*zh*(1 - 5*zh**2))
        gy=3.731763325901154E-1_r8*(3*ir*yh*zh*(1 - 5*zh**2))
        gz=3.731763325901154E-1_r8*(-3*ir*(1 - 6*zh**2 + 5*zh**4))
        gxx=3.731763325901154E-1_r8*(3*ir**2*zh*(1 - 5*zh**2 + xh**2*(-3 + 25*zh**2)))
        gxy=3.731763325901154E-1_r8*(3*ir**2*xh*yh*zh*(-3 + 25*zh**2))
        gxz=3.731763325901154E-1_r8*(3*ir**2*xh*(1 - 18*zh**2 + 25*zh**4))
        gyy=3.731763325901154E-1_r8*(3*ir**2*zh*(1 - 5*zh**2 + yh**2*(-3 + 25*zh**2)))
        gyz=3.731763325901154E-1_r8*(3*ir**2*yh*(1 - 18*zh**2 + 25*zh**4))
        gzz=3.731763325901154E-1_r8*(3*ir**2*zh*(13 - 38*zh**2 + 25*zh**4))
    case(14)
        ! l=3, m=1, iang=14 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[21/(2*Pi)]*xh*(-1 + 5*zh^2))/4
        ! dYlm/dx = (Sqrt[21/(2*Pi)]*(yh^2 - 4*zh^2 - 3*xh^2*(-1 + xh^2 + yh^2 - 4*zh^2)))/(4*r)
        ! dYlm/dy = (Sqrt[21/(2*Pi)]*xh*yh*(2 - 3*xh^2 - 3*yh^2 + 12*zh^2))/(4*r)
        ! dYlm/dz = -(Sqrt[21/(2*Pi)]*xh*zh*(8 + 3*xh^2 + 3*yh^2 - 12*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (3*Sqrt[21/(2*Pi)]*xh*(2 + 5*xh^4 - 3*yh^2 + 12*zh^2 + xh^2*(-7 + 5*yh^2 - 20*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (Sqrt[21/(2*Pi)]*yh*(2 + 15*xh^4 - 3*yh^2 + 12*zh^2 + 15*xh^2*(-1 + yh^2 - 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[21/(2*Pi)]*zh*(8 - 15*xh^4 + 3*yh^2 - 12*zh^2 - 15*xh^2*(1 + yh^2 - 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dydx = (Sqrt[21/(2*Pi)]*yh*(2 + 15*xh^4 - 3*yh^2 + 12*zh^2 + 15*xh^2*(-1 + yh^2 - 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dydy = (Sqrt[21/(2*Pi)]*xh*(2 + 15*yh^4 + 3*xh^2*(-1 + 5*yh^2) + 12*zh^2 - 15*yh^2*(1 + 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[21/(2*Pi)]*xh*yh*zh*(6 + 5*xh^2 + 5*yh^2 - 20*zh^2))/(4*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[21/(2*Pi)]*zh*(8 - 15*xh^4 + 3*yh^2 - 12*zh^2 - 15*xh^2*(1 + yh^2 - 4*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[21/(2*Pi)]*xh*yh*zh*(6 + 5*xh^2 + 5*yh^2 - 20*zh^2))/(4*r^2)
        ! d^2Ylm/dzdz = (Sqrt[21/(2*Pi)]*xh*(-8 + 60*zh^2 - 60*zh^4 + 3*xh^2*(-1 + 5*zh^2) + 3*yh^2*(-1 + 5*zh^2)))/(4*r^2)
        y=-4.570457994644657E-1_r8*(xh*(-1 + 5*zh**2))
        gx=-4.570457994644657E-1_r8*(-(ir*(yh**2 - 4*zh**2 - 3*xh**2*(-1 + xh**2 + yh**2 - 4*zh**2))))
        gy=-4.570457994644657E-1_r8*(ir*xh*yh*(-2 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gz=-4.570457994644657E-1_r8*(ir*xh*zh*(8 + 3*xh**2 + 3*yh**2 - 12*zh**2))
        gxx=-4.570457994644657E-1_r8*(-3*ir**2*xh*(2 + 5*xh**4 - 3*yh**2 + 12*zh**2 + xh**2*(-7 + 5*yh**2 - 20*zh**2)))
        gxy=-4.570457994644657E-1_r8*(-(ir**2*yh*(2 + 15*xh**4 - 3*yh**2 + 12*zh**2 + 15*xh**2*(-1 + yh**2 - 4*zh**2))))
        gxz=-4.570457994644657E-1_r8*(ir**2*zh*(8 - 15*xh**4 + 3*yh**2 - 12*zh**2 - 15*xh**2*(1 + yh**2 - 4*zh**2)))
        gyy=-4.570457994644657E-1_r8*(-(ir**2*xh*(2 + 15*yh**4 + 3*xh**2*(-1 + 5*yh**2) + 12*zh**2 - 15*yh**2*(1 + 4*zh**2))))
        gyz=-4.570457994644657E-1_r8*(-3*ir**2*xh*yh*zh*(6 + 5*xh**2 + 5*yh**2 - 20*zh**2))
        gzz=-4.570457994644657E-1_r8*(-(ir**2*xh*(-8 + 60*zh**2 - 60*zh**4 + 3*xh**2*(-1 + 5*zh**2) + 3*yh**2*(-1 + 5*zh**2))))
    case(15)
        ! l=3, m=2, iang=15 (=l^2+l+1+m)
        ! Ylm= (Sqrt[105/Pi]*(xh^2 - yh^2)*zh)/4
        ! dYlm/dx = (Sqrt[105/Pi]*xh*(2 - 3*xh^2 + 3*yh^2)*zh)/(4*r)
        ! dYlm/dy = (Sqrt[105/Pi]*yh*(-2 - 3*xh^2 + 3*yh^2)*zh)/(4*r)
        ! dYlm/dz = -(Sqrt[105/Pi]*(xh - yh)*(xh + yh)*(-1 + 3*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (Sqrt[105/Pi]*(2 + 15*xh^4 + 3*yh^2 - 15*xh^2*(1 + yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dxdy = (15*Sqrt[105/Pi]*xh*yh*(xh^2 - yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dxdz = (Sqrt[105/Pi]*xh*(2 - 6*zh^2 + yh^2*(3 - 15*zh^2) + 3*xh^2*(-1 + 5*zh^2)))/(4*r^2)
        ! d^2Ylm/dydx = (15*Sqrt[105/Pi]*xh*yh*(xh^2 - yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dydy = (Sqrt[105/Pi]*(-2 + 15*yh^2 - 15*yh^4 + 3*xh^2*(-1 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dydz = -(Sqrt[105/Pi]*yh*(2 - 6*zh^2 + xh^2*(3 - 15*zh^2) + 3*yh^2*(-1 + 5*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdx = (Sqrt[105/Pi]*xh*(2 - 6*zh^2 + yh^2*(3 - 15*zh^2) + 3*xh^2*(-1 + 5*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdy = -(Sqrt[105/Pi]*yh*(2 - 6*zh^2 + xh^2*(3 - 15*zh^2) + 3*yh^2*(-1 + 5*zh^2)))/(4*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[105/Pi]*(xh^2 - yh^2)*zh*(-3 + 5*zh^2))/(4*r^2)
        y=1.445305721320277E0_r8*((xh**2 - yh**2)*zh)
        gx=1.445305721320277E0_r8*(ir*xh*(2 - 3*xh**2 + 3*yh**2)*zh)
        gy=1.445305721320277E0_r8*(ir*yh*(-2 - 3*xh**2 + 3*yh**2)*zh)
        gz=1.445305721320277E0_r8*(-(ir*(xh - yh)*(xh + yh)*(-1 + 3*zh**2)))
        gxx=1.445305721320277E0_r8*(ir**2*(2 + 15*xh**4 + 3*yh**2 - 15*xh**2*(1 + yh**2))*zh)
        gxy=1.445305721320277E0_r8*(15*ir**2*xh*yh*(xh**2 - yh**2)*zh)
        gxz=1.445305721320277E0_r8*(ir**2*xh*(2 - 6*zh**2 + yh**2*(3 - 15*zh**2) + 3*xh**2*(-1 + 5*zh**2)))
        gyy=1.445305721320277E0_r8*(ir**2*(-2 + 15*yh**2 - 15*yh**4 + 3*xh**2*(-1 + 5*yh**2))*zh)
        gyz=1.445305721320277E0_r8*(-(ir**2*yh*(2 - 6*zh**2 + xh**2*(3 - 15*zh**2) + 3*yh**2*(-1 + 5*zh**2))))
        gzz=1.445305721320277E0_r8*(3*ir**2*(xh**2 - yh**2)*zh*(-3 + 5*zh**2))
    case(16)
        ! l=3, m=3, iang=16 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*(xh^4 + yh^2 - xh^2*(1 + 3*yh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*xh*yh*(2 + xh^2 - 3*yh^2))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh)/(4*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[35/(2*Pi)]*xh*(2 + 5*xh^4 + 9*yh^2 - xh^2*(7 + 15*yh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[35/(2*Pi)]*yh*(2 - 5*xh^4 - 3*yh^2 + 3*xh^2*(-1 + 5*yh^2)))/(4*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[35/(2*Pi)]*(5*xh^4 + 3*yh^2 - 3*xh^2*(1 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[35/(2*Pi)]*yh*(2 - 5*xh^4 - 3*yh^2 + 3*xh^2*(-1 + 5*yh^2)))/(4*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[35/(2*Pi)]*xh*(-2 + 15*yh^2 - 15*yh^4 + xh^2*(-1 + 5*yh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[35/(2*Pi)]*xh*yh*(6 + 5*xh^2 - 15*yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dzdx = (-3*Sqrt[35/(2*Pi)]*(5*xh^4 + 3*yh^2 - 3*xh^2*(1 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[35/(2*Pi)]*xh*yh*(6 + 5*xh^2 - 15*yh^2)*zh)/(4*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 5*zh^2))/(4*r^2)
        y=-5.900435899266435E-1_r8*(xh**3 - 3*xh*yh**2)
        gx=-5.900435899266435E-1_r8*(-3*ir*(xh**4 + yh**2 - xh**2*(1 + 3*yh**2)))
        gy=-5.900435899266435E-1_r8*(-3*ir*xh*yh*(2 + xh**2 - 3*yh**2))
        gz=-5.900435899266435E-1_r8*(-3*ir*xh*(xh**2 - 3*yh**2)*zh)
        gxx=-5.900435899266435E-1_r8*(3*ir**2*xh*(2 + 5*xh**4 + 9*yh**2 - xh**2*(7 + 15*yh**2)))
        gxy=-5.900435899266435E-1_r8*(-3*ir**2*yh*(2 - 5*xh**4 - 3*yh**2 + 3*xh**2*(-1 + 5*yh**2)))
        gxz=-5.900435899266435E-1_r8*(3*ir**2*(5*xh**4 + 3*yh**2 - 3*xh**2*(1 + 5*yh**2))*zh)
        gyy=-5.900435899266435E-1_r8*(3*ir**2*xh*(-2 + 15*yh**2 - 15*yh**4 + xh**2*(-1 + 5*yh**2)))
        gyz=-5.900435899266435E-1_r8*(3*ir**2*xh*yh*(6 + 5*xh**2 - 15*yh**2)*zh)
        gzz=-5.900435899266435E-1_r8*(3*ir**2*xh*(xh**2 - 3*yh**2)*(-1 + 5*zh**2))
    case(17)
        ! l=4, m=-4, iang=17 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*xh*yh*(xh^2 - yh^2))/4
        ! dYlm/dx = (3*Sqrt[35/Pi]*yh*(-4*xh^4 - yh^2 + xh^2*(3 + 4*yh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[35/Pi]*xh*(xh^2 - (3 + 4*xh^2)*yh^2 + 4*yh^4))/(4*r)
        ! dYlm/dz = (3*Sqrt[35/Pi]*xh*yh*(-xh^2 + yh^2)*zh)/r
        ! d^2Ylm/dxdx = (3*Sqrt[35/Pi]*xh*yh*(3 + 12*xh^4 + 6*yh^2 - 2*xh^2*(7 + 6*yh^2)))/(2*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[35/Pi]*(yh^2*(-3 + 4*yh^2) + 4*xh^4*(-1 + 6*yh^2) + xh^2*(3 - 24*yh^4)))/(4*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[35/Pi]*yh*(6*xh^4 + yh^2 - 3*xh^2*(1 + 2*yh^2))*zh)/r^2
        ! d^2Ylm/dydx = (3*Sqrt[35/Pi]*(yh^2*(-3 + 4*yh^2) + 4*xh^4*(-1 + 6*yh^2) + xh^2*(3 - 24*yh^4)))/(4*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[35/Pi]*xh*yh*(-3 + 14*yh^2 - 12*yh^4 + 6*xh^2*(-1 + 2*yh^2)))/(2*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[35/Pi]*xh*(3*yh^2 - 6*yh^4 + xh^2*(-1 + 6*yh^2))*zh)/r^2
        ! d^2Ylm/dzdx = (3*Sqrt[35/Pi]*yh*(6*xh^4 + yh^2 - 3*xh^2*(1 + 2*yh^2))*zh)/r^2
        ! d^2Ylm/dzdy = (3*Sqrt[35/Pi]*xh*(3*yh^2 - 6*yh^4 + xh^2*(-1 + 6*yh^2))*zh)/r^2
        ! d^2Ylm/dzdz = (3*Sqrt[35/Pi]*xh*yh*(xh^2 - yh^2)*(-1 + 6*zh^2))/r^2
        y=2.503342941796705E0_r8*(xh*yh*(xh**2 - yh**2))
        gx=2.503342941796705E0_r8*(ir*yh*(-4*xh**4 - yh**2 + xh**2*(3 + 4*yh**2)))
        gy=2.503342941796705E0_r8*(ir*xh*(xh**2 - (3 + 4*xh**2)*yh**2 + 4*yh**4))
        gz=2.503342941796705E0_r8*(4*ir*xh*yh*(-xh**2 + yh**2)*zh)
        gxx=2.503342941796705E0_r8*(2*ir**2*xh*yh*(3 + 12*xh**4 + 6*yh**2 - 2*xh**2*(7 + 6*yh**2)))
        gxy=2.503342941796705E0_r8*(ir**2*(yh**2*(-3 + 4*yh**2) + 4*xh**4*(-1 + 6*yh**2) + xh**2*(3 - 24*yh**4)))
        gxz=2.503342941796705E0_r8*(4*ir**2*yh*(6*xh**4 + yh**2 - 3*xh**2*(1 + 2*yh**2))*zh)
        gyy=2.503342941796705E0_r8*(2*ir**2*xh*yh*(-3 + 14*yh**2 - 12*yh**4 + 6*xh**2*(-1 + 2*yh**2)))
        gyz=2.503342941796705E0_r8*(4*ir**2*xh*(3*yh**2 - 6*yh**4 + xh**2*(-1 + 6*yh**2))*zh)
        gzz=2.503342941796705E0_r8*(4*ir**2*xh*yh*(xh**2 - yh**2)*(-1 + 6*zh**2))
    case(18)
        ! l=4, m=-3, iang=18 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 - 6*xh^2 + 2*yh^2)*zh)/(2*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*(3*xh^2 - 3*(1 + 4*xh^2)*yh^2 + 4*yh^4)*zh)/(4*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 4*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[35/(2*Pi)]*yh*(-3 - 36*xh^4 - 2*yh^2 + 6*xh^2*(5 + 2*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dxdy = (9*Sqrt[35/(2*Pi)]*xh*(1 - 2*yh^2 - 4*yh^4 + 2*xh^2*(-1 + 6*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 - 12*zh^2 + yh^2*(2 - 12*zh^2) + 6*xh^2*(-1 + 6*zh^2)))/(2*r^2)
        ! d^2Ylm/dydx = (9*Sqrt[35/(2*Pi)]*xh*(1 - 2*yh^2 - 4*yh^4 + 2*xh^2*(-1 + 6*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[35/(2*Pi)]*yh*(3 - 14*yh^2 + 12*yh^4 + xh^2*(18 - 36*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[35/(2*Pi)]*(yh^2*(-3 + 12*zh^2 + yh^2*(4 - 24*zh^2)) + 3*xh^2*(1 - 4*zh^2 + 4*yh^2*(-1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 - 12*zh^2 + yh^2*(2 - 12*zh^2) + 6*xh^2*(-1 + 6*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[35/(2*Pi)]*(yh^2*(-3 + 12*zh^2 + yh^2*(4 - 24*zh^2)) + 3*xh^2*(1 - 4*zh^2 + 4*yh^2*(-1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdz = (-9*Sqrt[35/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(-1 + 2*zh^2))/r^2
        y=-1.770130769779931E0_r8*(yh*(-3*xh**2 + yh**2)*zh)
        gx=-1.770130769779931E0_r8*(2*ir*xh*yh*(-3 + 6*xh**2 - 2*yh**2)*zh)
        gy=-1.770130769779931E0_r8*(-(ir*(3*xh**2 - 3*(1 + 4*xh**2)*yh**2 + 4*yh**4)*zh))
        gz=-1.770130769779931E0_r8*(-(ir*yh*(-3*xh**2 + yh**2)*(-1 + 4*zh**2)))
        gxx=-1.770130769779931E0_r8*(2*ir**2*yh*(-3 - 36*xh**4 - 2*yh**2 + 6*xh**2*(5 + 2*yh**2))*zh)
        gxy=-1.770130769779931E0_r8*(-6*ir**2*xh*(1 - 2*yh**2 - 4*yh**4 + 2*xh**2*(-1 + 6*yh**2))*zh)
        gxz=-1.770130769779931E0_r8*(-2*ir**2*xh*yh*(3 - 12*zh**2 + yh**2*(2 - 12*zh**2) + 6*xh**2*(-1 + 6*zh**2)))
        gyy=-1.770130769779931E0_r8*(2*ir**2*yh*(3 - 14*yh**2 + 12*yh**4 + xh**2*(18 - 36*yh**2))*zh)
        gyz=-1.770130769779931E0_r8*(-(ir**2*(yh**2*(-3 + 12*zh**2 + yh**2*(4 - 24*zh**2)) + 3*xh**2*(1 - 4*zh**2 + 4*yh**2*(-1 + 6*zh**2)))))
        gzz=-1.770130769779931E0_r8*(12*ir**2*yh*(-3*xh**2 + yh**2)*zh*(-1 + 2*zh**2))
    case(19)
        ! l=4, m=-2, iang=19 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*xh*yh*(-1 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/Pi]*yh*(4*xh^4 - yh^2 + 6*zh^2 + xh^2*(-3 + 4*yh^2 - 24*zh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[5/Pi]*xh*(4*yh^4 + xh^2*(-1 + 4*yh^2) + 6*zh^2 - 3*yh^2*(1 + 8*zh^2)))/(4*r)
        ! dYlm/dz = (3*Sqrt[5/Pi]*xh*yh*zh*(3 + xh^2 + yh^2 - 6*zh^2))/r
        ! d^2Ylm/dxdx = (-3*Sqrt[5/Pi]*xh*yh*(3 + 12*xh^4 - 6*yh^2 + 36*zh^2 + 2*xh^2*(-7 + 6*yh^2 - 36*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[5/Pi]*(-4*yh^4 + 4*xh^4*(-1 + 6*yh^2) - 6*zh^2 + 3*yh^2*(1 + 8*zh^2) + 3*xh^2*(1 + 8*yh^4 + 8*zh^2 - 8*yh^2*(1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[5/Pi]*yh*zh*(-3 + 6*xh^4 - yh^2 + 6*zh^2 + xh^2*(9 + 6*yh^2 - 36*zh^2)))/r^2
        ! d^2Ylm/dydx = (-3*Sqrt[5/Pi]*(-4*yh^4 + 4*xh^4*(-1 + 6*yh^2) - 6*zh^2 + 3*yh^2*(1 + 8*zh^2) + 3*xh^2*(1 + 8*yh^4 + 8*zh^2 - 8*yh^2*(1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[5/Pi]*xh*yh*(3 + 12*yh^4 + 6*xh^2*(-1 + 2*yh^2) + 36*zh^2 - 2*yh^2*(7 + 36*zh^2)))/(2*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[5/Pi]*xh*zh*(-3 + 6*yh^4 + xh^2*(-1 + 6*yh^2) + 6*zh^2 + yh^2*(9 - 36*zh^2)))/r^2
        ! d^2Ylm/dzdx = (-3*Sqrt[5/Pi]*yh*zh*(-3 + 6*xh^4 - yh^2 + 6*zh^2 + xh^2*(9 + 6*yh^2 - 36*zh^2)))/r^2
        ! d^2Ylm/dzdy = (-3*Sqrt[5/Pi]*xh*zh*(-3 + 6*yh^4 + xh^2*(-1 + 6*yh^2) + 6*zh^2 + yh^2*(9 - 36*zh^2)))/r^2
        ! d^2Ylm/dzdz = (3*Sqrt[5/Pi]*xh*yh*(3 - 30*zh^2 + 36*zh^4 + xh^2*(1 - 6*zh^2) + yh^2*(1 - 6*zh^2)))/r^2
        y=9.461746957575600E-1_r8*(xh*yh*(-1 + 7*zh**2))
        gx=9.461746957575600E-1_r8*(ir*yh*(4*xh**4 - yh**2 + 6*zh**2 + xh**2*(-3 + 4*yh**2 - 24*zh**2)))
        gy=9.461746957575600E-1_r8*(ir*xh*(4*yh**4 + xh**2*(-1 + 4*yh**2) + 6*zh**2 - 3*yh**2*(1 + 8*zh**2)))
        gz=9.461746957575600E-1_r8*(4*ir*xh*yh*zh*(3 + xh**2 + yh**2 - 6*zh**2))
        gxx=9.461746957575600E-1_r8*(-2*ir**2*xh*yh*(3 + 12*xh**4 - 6*yh**2 + 36*zh**2 + 2*xh**2*(-7 + 6*yh**2 - 36*zh**2)))
        gxy=9.461746957575600E-1_r8*(-(ir**2*(-4*yh**4 + 4*xh**4*(-1 + 6*yh**2) - 6*zh**2 + 3*yh**2*(1 + 8*zh**2) + 3*xh**2*(1 + 8*yh**4 + 8*zh**2 - 8*yh**2*(1 + 6*zh**2)))))
        gxz=9.461746957575600E-1_r8*(-4*ir**2*yh*zh*(-3 + 6*xh**4 - yh**2 + 6*zh**2 + xh**2*(9 + 6*yh**2 - 36*zh**2)))
        gyy=9.461746957575600E-1_r8*(-2*ir**2*xh*yh*(3 + 12*yh**4 + 6*xh**2*(-1 + 2*yh**2) + 36*zh**2 - 2*yh**2*(7 + 36*zh**2)))
        gyz=9.461746957575600E-1_r8*(-4*ir**2*xh*zh*(-3 + 6*yh**4 + xh**2*(-1 + 6*yh**2) + 6*zh**2 + yh**2*(9 - 36*zh**2)))
        gzz=9.461746957575600E-1_r8*(4*ir**2*xh*yh*(3 - 30*zh**2 + 36*zh**4 + xh**2*(1 - 6*zh**2) + yh**2*(1 - 6*zh**2)))
    case(20)
        ! l=4, m=-1, iang=20 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/(2*Pi)]*yh*zh*(-3 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/(2*Pi)]*xh*yh*zh*(3 - 14*zh^2))/(2*r)
        ! dYlm/dy = (3*Sqrt[5/(2*Pi)]*zh*(-3 + 7*zh^2 + yh^2*(6 - 28*zh^2)))/(4*r)
        ! dYlm/dz = (-3*Sqrt[5/(2*Pi)]*yh*(3 - 27*zh^2 + 28*zh^4))/(4*r)
        ! d^2Ylm/dxdx = (3*Sqrt[5/(2*Pi)]*yh*zh*(3 - 14*zh^2 + 12*xh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[5/(2*Pi)]*xh*zh*(3 - 14*zh^2 + 12*yh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdz = (9*Sqrt[5/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 28*zh^4))/(2*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[5/(2*Pi)]*xh*zh*(3 - 14*zh^2 + 12*yh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dydy = (9*Sqrt[5/(2*Pi)]*yh*zh*(3 - 14*zh^2 + 4*yh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[5/(2*Pi)]*(-3 + 27*zh^2 - 28*zh^4 + 6*yh^2*(1 - 18*zh^2 + 28*zh^4)))/(4*r^2)
        ! d^2Ylm/dzdx = (9*Sqrt[5/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 28*zh^4))/(2*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[5/(2*Pi)]*(-3 + 27*zh^2 - 28*zh^4 + 6*yh^2*(1 - 18*zh^2 + 28*zh^4)))/(4*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[5/(2*Pi)]*yh*zh*(15 - 55*zh^2 + 42*zh^4))/r^2
        y=6.690465435572892E-1_r8*(yh*zh*(-3 + 7*zh**2))
        gx=6.690465435572892E-1_r8*(2*ir*xh*yh*zh*(3 - 14*zh**2))
        gy=6.690465435572892E-1_r8*(ir*zh*(-3 + 7*zh**2 + yh**2*(6 - 28*zh**2)))
        gz=6.690465435572892E-1_r8*(-(ir*yh*(3 - 27*zh**2 + 28*zh**4)))
        gxx=6.690465435572892E-1_r8*(2*ir**2*yh*zh*(3 - 14*zh**2 + 12*xh**2*(-1 + 7*zh**2)))
        gxy=6.690465435572892E-1_r8*(2*ir**2*xh*zh*(3 - 14*zh**2 + 12*yh**2*(-1 + 7*zh**2)))
        gxz=6.690465435572892E-1_r8*(6*ir**2*xh*yh*(1 - 18*zh**2 + 28*zh**4))
        gyy=6.690465435572892E-1_r8*(6*ir**2*yh*zh*(3 - 14*zh**2 + 4*yh**2*(-1 + 7*zh**2)))
        gyz=6.690465435572892E-1_r8*(ir**2*(-3 + 27*zh**2 - 28*zh**4 + 6*yh**2*(1 - 18*zh**2 + 28*zh**4)))
        gzz=6.690465435572892E-1_r8*(4*ir**2*yh*zh*(15 - 55*zh**2 + 42*zh**4))
    case(21)
        ! l=4, m=0, iang=21 (=l^2+l+1+m)
        ! Ylm= (3*(3 - 30*zh^2 + 35*zh^4))/(16*Sqrt[Pi])
        ! dYlm/dx = (15*xh*zh^2*(3 - 7*zh^2))/(4*Sqrt[Pi]*r)
        ! dYlm/dy = (15*yh*zh^2*(3 - 7*zh^2))/(4*Sqrt[Pi]*r)
        ! dYlm/dz = (-15*zh*(3 - 10*zh^2 + 7*zh^4))/(4*Sqrt[Pi]*r)
        ! d^2Ylm/dxdx = (15*zh^2*(3 - 7*zh^2 + 6*xh^2*(-2 + 7*zh^2)))/(4*Sqrt[Pi]*r^2)
        ! d^2Ylm/dxdy = (45*xh*yh*zh^2*(-2 + 7*zh^2))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dxdz = (15*xh*zh*(3 - 20*zh^2 + 21*zh^4))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dydx = (45*xh*yh*zh^2*(-2 + 7*zh^2))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dydy = (15*zh^2*(3 - 7*zh^2 + 6*yh^2*(-2 + 7*zh^2)))/(4*Sqrt[Pi]*r^2)
        ! d^2Ylm/dydz = (15*yh*zh*(3 - 20*zh^2 + 21*zh^4))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dzdx = (15*xh*zh*(3 - 20*zh^2 + 21*zh^4))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dzdy = (15*yh*zh*(3 - 20*zh^2 + 21*zh^4))/(2*Sqrt[Pi]*r^2)
        ! d^2Ylm/dzdz = (45*(-1 + 12*zh^2 - 25*zh^4 + 14*zh^6))/(4*Sqrt[Pi]*r^2)
        y=1.057855469152043E-1_r8*(3 - 30*zh**2 + 35*zh**4)
        gx=1.057855469152043E-1_r8*(20*ir*xh*zh**2*(3 - 7*zh**2))
        gy=1.057855469152043E-1_r8*(20*ir*yh*zh**2*(3 - 7*zh**2))
        gz=1.057855469152043E-1_r8*(-20*ir*zh*(3 - 10*zh**2 + 7*zh**4))
        gxx=1.057855469152043E-1_r8*(20*ir**2*zh**2*(3 - 7*zh**2 + 6*xh**2*(-2 + 7*zh**2)))
        gxy=1.057855469152043E-1_r8*(120*ir**2*xh*yh*zh**2*(-2 + 7*zh**2))
        gxz=1.057855469152043E-1_r8*(40*ir**2*xh*zh*(3 - 20*zh**2 + 21*zh**4))
        gyy=1.057855469152043E-1_r8*(20*ir**2*zh**2*(3 - 7*zh**2 + 6*yh**2*(-2 + 7*zh**2)))
        gyz=1.057855469152043E-1_r8*(40*ir**2*yh*zh*(3 - 20*zh**2 + 21*zh**4))
        gzz=1.057855469152043E-1_r8*(60*ir**2*(-1 + 12*zh**2 - 25*zh**4 + 14*zh**6))
    case(22)
        ! l=4, m=1, iang=22 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[5/(2*Pi)]*xh*zh*(-3 + 7*zh^2))/4
        ! dYlm/dx = (3*Sqrt[5/(2*Pi)]*zh*(3 - 7*zh^2 + xh^2*(-6 + 28*zh^2)))/(4*r)
        ! dYlm/dy = (3*Sqrt[5/(2*Pi)]*xh*yh*zh*(-3 + 14*zh^2))/(2*r)
        ! dYlm/dz = (3*Sqrt[5/(2*Pi)]*xh*(3 - 27*zh^2 + 28*zh^4))/(4*r)
        ! d^2Ylm/dxdx = (-9*Sqrt[5/(2*Pi)]*xh*zh*(3 - 14*zh^2 + 4*xh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[5/(2*Pi)]*yh*zh*(3 - 14*zh^2 + 12*xh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[5/(2*Pi)]*(3 - 27*zh^2 + 28*zh^4 - 6*xh^2*(1 - 18*zh^2 + 28*zh^4)))/(4*r^2)
        ! d^2Ylm/dydx = (-3*Sqrt[5/(2*Pi)]*yh*zh*(3 - 14*zh^2 + 12*xh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[5/(2*Pi)]*xh*zh*(3 - 14*zh^2 + 12*yh^2*(-1 + 7*zh^2)))/(2*r^2)
        ! d^2Ylm/dydz = (-9*Sqrt[5/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 28*zh^4))/(2*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[5/(2*Pi)]*(3 - 27*zh^2 + 28*zh^4 - 6*xh^2*(1 - 18*zh^2 + 28*zh^4)))/(4*r^2)
        ! d^2Ylm/dzdy = (-9*Sqrt[5/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 28*zh^4))/(2*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[5/(2*Pi)]*xh*zh*(15 - 55*zh^2 + 42*zh^4))/r^2
        y=-6.690465435572892E-1_r8*(xh*zh*(-3 + 7*zh**2))
        gx=-6.690465435572892E-1_r8*(-(ir*zh*(3 - 7*zh**2 + xh**2*(-6 + 28*zh**2))))
        gy=-6.690465435572892E-1_r8*(-2*ir*xh*yh*zh*(-3 + 14*zh**2))
        gz=-6.690465435572892E-1_r8*(-(ir*xh*(3 - 27*zh**2 + 28*zh**4)))
        gxx=-6.690465435572892E-1_r8*(6*ir**2*xh*zh*(3 - 14*zh**2 + 4*xh**2*(-1 + 7*zh**2)))
        gxy=-6.690465435572892E-1_r8*(2*ir**2*yh*zh*(3 - 14*zh**2 + 12*xh**2*(-1 + 7*zh**2)))
        gxz=-6.690465435572892E-1_r8*(-(ir**2*(3 - 27*zh**2 + 28*zh**4 - 6*xh**2*(1 - 18*zh**2 + 28*zh**4))))
        gyy=-6.690465435572892E-1_r8*(2*ir**2*xh*zh*(3 - 14*zh**2 + 12*yh**2*(-1 + 7*zh**2)))
        gyz=-6.690465435572892E-1_r8*(6*ir**2*xh*yh*(1 - 18*zh**2 + 28*zh**4))
        gzz=-6.690465435572892E-1_r8*(4*ir**2*xh*zh*(15 - 55*zh**2 + 42*zh**4))
    case(23)
        ! l=4, m=2, iang=23 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[5/Pi]*(xh^2 - yh^2)*(-1 + 7*zh^2))/8
        ! dYlm/dx = (3*Sqrt[5/Pi]*xh*(-xh^2 + xh^4 - yh^4 + 3*(1 - 2*xh^2 + 2*yh^2)*zh^2))/(2*r)
        ! dYlm/dy = (3*Sqrt[5/Pi]*yh*(xh^4 + yh^2 - yh^4 - 3*(1 + 2*xh^2 - 2*yh^2)*zh^2))/(2*r)
        ! dYlm/dz = (3*Sqrt[5/Pi]*(xh - yh)*(xh + yh)*zh*(3 + xh^2 + yh^2 - 6*zh^2))/(2*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[5/Pi]*(6*xh^6 + yh^4 - 3*zh^2 - 6*yh^2*zh^2 - 9*xh^4*(1 + 4*zh^2) + xh^2*(3 - 6*yh^4 + 30*zh^2 + 36*yh^2*zh^2)))/(2*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[5/Pi]*xh*yh*(xh^2 - yh^2)*(-2 + 3*xh^2 + 3*yh^2 - 18*zh^2))/r^2
        ! d^2Ylm/dxdz = (-3*Sqrt[5/Pi]*xh*zh*(3*xh^4 + xh^2*(4 - 18*zh^2) - 3*(1 + yh^4 - 2*zh^2 + yh^2*(2 - 6*zh^2))))/r^2
        ! d^2Ylm/dydx = (-3*Sqrt[5/Pi]*xh*yh*(xh^2 - yh^2)*(-2 + 3*xh^2 + 3*yh^2 - 18*zh^2))/r^2
        ! d^2Ylm/dydy = (3*Sqrt[5/Pi]*(6*yh^6 + xh^4*(1 - 6*yh^2) - 3*zh^2 + 6*xh^2*(-1 + 6*yh^2)*zh^2 - 9*yh^4*(1 + 4*zh^2) + yh^2*(3 + 30*zh^2)))/(2*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[5/Pi]*yh*zh*(-3 - 3*xh^4 + 3*yh^4 + 6*zh^2 + yh^2*(4 - 18*zh^2) + 6*xh^2*(-1 + 3*zh^2)))/r^2
        ! d^2Ylm/dzdx = (-3*Sqrt[5/Pi]*xh*zh*(3*xh^4 + xh^2*(4 - 18*zh^2) - 3*(1 + yh^4 - 2*zh^2 + yh^2*(2 - 6*zh^2))))/r^2
        ! d^2Ylm/dzdy = (3*Sqrt[5/Pi]*yh*zh*(-3 - 3*xh^4 + 3*yh^4 + 6*zh^2 + yh^2*(4 - 18*zh^2) + 6*xh^2*(-1 + 3*zh^2)))/r^2
        ! d^2Ylm/dzdz = (-3*Sqrt[5/Pi]*(xh^2 - yh^2)*(-3 + 30*zh^2 - 36*zh^4 + xh^2*(-1 + 6*zh^2) + yh^2*(-1 + 6*zh^2)))/(2*r^2)
        y=4.730873478787800E-1_r8*((xh**2 - yh**2)*(-1 + 7*zh**2))
        gx=4.730873478787800E-1_r8*(4*ir*xh*(-xh**2 + xh**4 - yh**4 + 3*(1 - 2*xh**2 + 2*yh**2)*zh**2))
        gy=4.730873478787800E-1_r8*(4*ir*yh*(xh**4 + yh**2 - yh**4 - 3*(1 + 2*xh**2 - 2*yh**2)*zh**2))
        gz=4.730873478787800E-1_r8*(4*ir*(xh - yh)*(xh + yh)*zh*(3 + xh**2 + yh**2 - 6*zh**2))
        gxx=4.730873478787800E-1_r8*(-4*ir**2*(6*xh**6 + yh**4 - 3*zh**2 - 6*yh**2*zh**2 - 9*xh**4*(1 + 4*zh**2) + xh**2*(3 - 6*yh**4 + 30*zh**2 + 36*yh**2*zh**2)))
        gxy=4.730873478787800E-1_r8*(-8*ir**2*xh*yh*(xh**2 - yh**2)*(-2 + 3*xh**2 + 3*yh**2 - 18*zh**2))
        gxz=4.730873478787800E-1_r8*(-8*ir**2*xh*zh*(3*xh**4 + xh**2*(4 - 18*zh**2) - 3*(1 + yh**4 - 2*zh**2 + yh**2*(2 - 6*zh**2))))
        gyy=4.730873478787800E-1_r8*(4*ir**2*(6*yh**6 + xh**4*(1 - 6*yh**2) - 3*zh**2 + 6*xh**2*(-1 + 6*yh**2)*zh**2 - 9*yh**4*(1 + 4*zh**2) + yh**2*(3 + 30*zh**2)))
        gyz=4.730873478787800E-1_r8*(8*ir**2*yh*zh*(-3 - 3*xh**4 + 3*yh**4 + 6*zh**2 + yh**2*(4 - 18*zh**2) + 6*xh**2*(-1 + 3*zh**2)))
        gzz=4.730873478787800E-1_r8*(-4*ir**2*(xh**2 - yh**2)*(-3 + 30*zh**2 - 36*zh**4 + xh**2*(-1 + 6*zh**2) + yh**2*(-1 + 6*zh**2)))
    case(24)
        ! l=4, m=3, iang=24 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[35/(2*Pi)]*(4*xh^4 + 3*yh^2 - 3*xh^2*(1 + 4*yh^2))*zh)/(4*r)
        ! dYlm/dy = (3*Sqrt[35/(2*Pi)]*xh*yh*(3 + 2*xh^2 - 6*yh^2)*zh)/(2*r)
        ! dYlm/dz = (3*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 4*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[35/(2*Pi)]*xh*(3 + 12*xh^4 + 18*yh^2 - 2*xh^2*(7 + 18*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dxdy = (9*Sqrt[35/(2*Pi)]*yh*(1 - 4*xh^4 - 2*yh^2 + 2*xh^2*(-1 + 6*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[35/(2*Pi)]*(xh^4*(4 - 24*zh^2) + 3*yh^2*(1 - 4*zh^2) + 3*xh^2*(-1 + 4*zh^2 + 4*yh^2*(-1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dydx = (9*Sqrt[35/(2*Pi)]*yh*(1 - 4*xh^4 - 2*yh^2 + 2*xh^2*(-1 + 6*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[35/(2*Pi)]*xh*(-3 + 30*yh^2 - 36*yh^4 + 2*xh^2*(-1 + 6*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[35/(2*Pi)]*xh*yh*(-3 + 12*zh^2 + yh^2*(6 - 36*zh^2) + 2*xh^2*(-1 + 6*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[35/(2*Pi)]*(xh^4*(4 - 24*zh^2) + 3*yh^2*(1 - 4*zh^2) + 3*xh^2*(-1 + 4*zh^2 + 4*yh^2*(-1 + 6*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[35/(2*Pi)]*xh*yh*(-3 + 12*zh^2 + yh^2*(6 - 36*zh^2) + 2*xh^2*(-1 + 6*zh^2)))/(2*r^2)
        ! d^2Ylm/dzdz = (-9*Sqrt[35/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(-1 + 2*zh^2))/r^2
        y=-1.770130769779931E0_r8*(xh*(xh**2 - 3*yh**2)*zh)
        gx=-1.770130769779931E0_r8*(-(ir*(4*xh**4 + 3*yh**2 - 3*xh**2*(1 + 4*yh**2))*zh))
        gy=-1.770130769779931E0_r8*(-2*ir*xh*yh*(3 + 2*xh**2 - 6*yh**2)*zh)
        gz=-1.770130769779931E0_r8*(-(ir*xh*(xh**2 - 3*yh**2)*(-1 + 4*zh**2)))
        gxx=-1.770130769779931E0_r8*(2*ir**2*xh*(3 + 12*xh**4 + 18*yh**2 - 2*xh**2*(7 + 18*yh**2))*zh)
        gxy=-1.770130769779931E0_r8*(-6*ir**2*yh*(1 - 4*xh**4 - 2*yh**2 + 2*xh**2*(-1 + 6*yh**2))*zh)
        gxz=-1.770130769779931E0_r8*(-(ir**2*(xh**4*(4 - 24*zh**2) + 3*yh**2*(1 - 4*zh**2) + 3*xh**2*(-1 + 4*zh**2 + 4*yh**2*(-1 + 6*zh**2)))))
        gyy=-1.770130769779931E0_r8*(2*ir**2*xh*(-3 + 30*yh**2 - 36*yh**4 + 2*xh**2*(-1 + 6*yh**2))*zh)
        gyz=-1.770130769779931E0_r8*(2*ir**2*xh*yh*(-3 + 12*zh**2 + yh**2*(6 - 36*zh**2) + 2*xh**2*(-1 + 6*zh**2)))
        gzz=-1.770130769779931E0_r8*(12*ir**2*xh*(xh**2 - 3*yh**2)*zh*(-1 + 2*zh**2))
    case(25)
        ! l=4, m=4, iang=25 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4))/16
        ! dYlm/dx = (-3*Sqrt[35/Pi]*xh*(xh^4 + 3*yh^2 + yh^4 - xh^2*(1 + 6*yh^2)))/(4*r)
        ! dYlm/dy = (-3*Sqrt[35/Pi]*yh*(xh^4 - yh^2 + yh^4 + xh^2*(3 - 6*yh^2)))/(4*r)
        ! dYlm/dz = (-3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh)/(4*r)
        ! d^2Ylm/dxdx = (3*Sqrt[35/Pi]*(6*xh^6 - yh^2*(3 + yh^2) - 9*xh^4*(1 + 4*yh^2) + xh^2*(3 + 30*yh^2 + 6*yh^4)))/(4*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[35/Pi]*xh*yh*(-3 + 3*xh^4 + 4*yh^2 + 3*yh^4 + xh^2*(4 - 18*yh^2)))/(2*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[35/Pi]*(3*xh^5 + 3*xh*yh^2*(2 + yh^2) - 2*xh^3*(1 + 9*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[35/Pi]*xh*yh*(-3 + 3*xh^4 + 4*yh^2 + 3*yh^4 + xh^2*(4 - 18*yh^2)))/(2*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[35/Pi]*(xh^4*(-1 + 6*yh^2) + xh^2*(-3 + 30*yh^2 - 36*yh^4) + 3*(yh^2 - 3*yh^4 + 2*yh^6)))/(4*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[35/Pi]*yh*(3*xh^4 + xh^2*(6 - 18*yh^2) + yh^2*(-2 + 3*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[35/Pi]*(3*xh^5 + 3*xh*yh^2*(2 + yh^2) - 2*xh^3*(1 + 9*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[35/Pi]*yh*(3*xh^4 + xh^2*(6 - 18*yh^2) + yh^2*(-2 + 3*yh^2))*zh)/(2*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[35/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 6*zh^2))/(4*r^2)
        y=6.258357354491761E-1_r8*(xh**4 - 6*xh**2*yh**2 + yh**4)
        gx=6.258357354491761E-1_r8*(-4*ir*xh*(xh**4 + 3*yh**2 + yh**4 - xh**2*(1 + 6*yh**2)))
        gy=6.258357354491761E-1_r8*(-4*ir*yh*(xh**4 - yh**2 + yh**4 + xh**2*(3 - 6*yh**2)))
        gz=6.258357354491761E-1_r8*(-4*ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*zh)
        gxx=6.258357354491761E-1_r8*(4*ir**2*(6*xh**6 - yh**2*(3 + yh**2) - 9*xh**4*(1 + 4*yh**2) + xh**2*(3 + 30*yh**2 + 6*yh**4)))
        gxy=6.258357354491761E-1_r8*(8*ir**2*xh*yh*(-3 + 3*xh**4 + 4*yh**2 + 3*yh**4 + xh**2*(4 - 18*yh**2)))
        gxz=6.258357354491761E-1_r8*(8*ir**2*(3*xh**5 + 3*xh*yh**2*(2 + yh**2) - 2*xh**3*(1 + 9*yh**2))*zh)
        gyy=6.258357354491761E-1_r8*(4*ir**2*(xh**4*(-1 + 6*yh**2) + xh**2*(-3 + 30*yh**2 - 36*yh**4) + 3*(yh**2 - 3*yh**4 + 2*yh**6)))
        gyz=6.258357354491761E-1_r8*(8*ir**2*yh*(3*xh**4 + xh**2*(6 - 18*yh**2) + yh**2*(-2 + 3*yh**2))*zh)
        gzz=6.258357354491761E-1_r8*(4*ir**2*(xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 6*zh**2))
    case(26)
        ! l=5, m=-5, iang=26 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4))/16
        ! dYlm/dx = (-15*Sqrt[77/(2*Pi)]*xh*yh*(5*xh^4 + 4*yh^2 + yh^4 - 2*xh^2*(2 + 5*yh^2)))/(16*r)
        ! dYlm/dy = (-15*Sqrt[77/(2*Pi)]*(-yh^4 + yh^6 + 2*xh^2*yh^2*(3 - 5*yh^2) + xh^4*(-1 + 5*yh^2)))/(16*r)
        ! dYlm/dz = (-15*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh)/(16*r)
        ! d^2Ylm/dxdx = (15*Sqrt[77/(2*Pi)]*yh*(35*xh^6 - yh^2*(4 + yh^2) - 5*xh^4*(9 + 14*yh^2) + xh^2*(12 + 50*yh^2 + 7*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (15*Sqrt[77/(2*Pi)]*(5*xh^5*(-1 + 7*yh^2) + xh^3*(4 + 10*yh^2 - 70*yh^4) + xh*yh^2*(-12 + 15*yh^2 + 7*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdz = (15*Sqrt[77/(2*Pi)]*yh*(35*xh^5 - 10*xh^3*(2 + 7*yh^2) + xh*yh^2*(20 + 7*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dydx = (15*Sqrt[77/(2*Pi)]*(5*xh^5*(-1 + 7*yh^2) + xh^3*(4 + 10*yh^2 - 70*yh^4) + xh*yh^2*(-12 + 15*yh^2 + 7*yh^4)))/(16*r^2)
        ! d^2Ylm/dydy = (15*Sqrt[77/(2*Pi)]*yh*(5*xh^4*(-3 + 7*yh^2) + yh^2*(4 - 11*yh^2 + 7*yh^4) - 2*xh^2*(6 - 35*yh^2 + 35*yh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (15*Sqrt[77/(2*Pi)]*(yh^4*(-5 + 7*yh^2) + 5*xh^4*(-1 + 7*yh^2) + xh^2*(30*yh^2 - 70*yh^4))*zh)/(16*r^2)
        ! d^2Ylm/dzdx = (15*Sqrt[77/(2*Pi)]*yh*(35*xh^5 - 10*xh^3*(2 + 7*yh^2) + xh*yh^2*(20 + 7*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dzdy = (15*Sqrt[77/(2*Pi)]*(yh^4*(-5 + 7*yh^2) + 5*xh^4*(-1 + 7*yh^2) + xh^2*(30*yh^2 - 70*yh^4))*zh)/(16*r^2)
        ! d^2Ylm/dzdz = (15*Sqrt[77/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*(-1 + 7*zh^2))/(16*r^2)
        y=6.563820568401701E-1_r8*(5*xh**4*yh - 10*xh**2*yh**3 + yh**5)
        gx=6.563820568401701E-1_r8*(-5*ir*xh*yh*(5*xh**4 + 4*yh**2 + yh**4 - 2*xh**2*(2 + 5*yh**2)))
        gy=6.563820568401701E-1_r8*(-5*ir*(-yh**4 + yh**6 + 2*xh**2*yh**2*(3 - 5*yh**2) + xh**4*(-1 + 5*yh**2)))
        gz=6.563820568401701E-1_r8*(-5*ir*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh)
        gxx=6.563820568401701E-1_r8*(5*ir**2*yh*(35*xh**6 - yh**2*(4 + yh**2) - 5*xh**4*(9 + 14*yh**2) + xh**2*(12 + 50*yh**2 + 7*yh**4)))
        gxy=6.563820568401701E-1_r8*(5*ir**2*(5*xh**5*(-1 + 7*yh**2) + xh**3*(4 + 10*yh**2 - 70*yh**4) + xh*yh**2*(-12 + 15*yh**2 + 7*yh**4)))
        gxz=6.563820568401701E-1_r8*(5*ir**2*yh*(35*xh**5 - 10*xh**3*(2 + 7*yh**2) + xh*yh**2*(20 + 7*yh**2))*zh)
        gyy=6.563820568401701E-1_r8*(5*ir**2*yh*(5*xh**4*(-3 + 7*yh**2) + yh**2*(4 - 11*yh**2 + 7*yh**4) - 2*xh**2*(6 - 35*yh**2 + 35*yh**4)))
        gyz=6.563820568401701E-1_r8*(5*ir**2*(yh**4*(-5 + 7*yh**2) + 5*xh**4*(-1 + 7*yh**2) + xh**2*(30*yh**2 - 70*yh**4))*zh)
        gzz=6.563820568401701E-1_r8*(5*ir**2*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*(-1 + 7*zh**2))
    case(27)
        ! l=5, m=-4, iang=27 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*xh*yh*(xh^2 - yh^2)*zh)/4
        ! dYlm/dx = (3*Sqrt[385/Pi]*yh*(-5*xh^4 - yh^2 + xh^2*(3 + 5*yh^2))*zh)/(4*r)
        ! dYlm/dy = (3*Sqrt[385/Pi]*xh*(xh^2 - (3 + 5*xh^2)*yh^2 + 5*yh^4)*zh)/(4*r)
        ! dYlm/dz = (-3*Sqrt[385/Pi]*xh*(xh - yh)*yh*(xh + yh)*(-1 + 5*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[385/Pi]*yh*(-35*xh^5 + 35*xh^3*(1 + yh^2) - 3*xh*(2 + 5*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[385/Pi]*(yh^2*(-3 + 5*yh^2) + 5*xh^4*(-1 + 7*yh^2) + xh^2*(3 - 35*yh^4))*zh)/(4*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[385/Pi]*yh*(xh^4*(5 - 35*zh^2) + yh^2*(1 - 5*zh^2) + xh^2*(-3 + 15*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[385/Pi]*(yh^2*(-3 + 5*yh^2) + 5*xh^4*(-1 + 7*yh^2) + xh^2*(3 - 35*yh^4))*zh)/(4*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[385/Pi]*xh*yh*(-6 + 35*yh^2 - 35*yh^4 + 5*xh^2*(-3 + 7*yh^2))*zh)/(4*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[385/Pi]*xh*(yh^2*(-3 + 15*zh^2 + yh^2*(5 - 35*zh^2)) + xh^2*(1 - 5*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdx = (-3*Sqrt[385/Pi]*yh*(xh^4*(5 - 35*zh^2) + yh^2*(1 - 5*zh^2) + xh^2*(-3 + 15*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[385/Pi]*xh*(yh^2*(-3 + 15*zh^2 + yh^2*(5 - 35*zh^2)) + xh^2*(1 - 5*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdz = (15*Sqrt[385/Pi]*xh*yh*(xh^2 - yh^2)*zh*(-3 + 7*zh^2))/(4*r^2)
        y=8.302649259524165E0_r8*(xh*yh*(xh**2 - yh**2)*zh)
        gx=8.302649259524165E0_r8*(ir*yh*(-5*xh**4 - yh**2 + xh**2*(3 + 5*yh**2))*zh)
        gy=8.302649259524165E0_r8*(ir*xh*(xh**2 - (3 + 5*xh**2)*yh**2 + 5*yh**4)*zh)
        gz=8.302649259524165E0_r8*(-(ir*xh*(xh - yh)*yh*(xh + yh)*(-1 + 5*zh**2)))
        gxx=8.302649259524165E0_r8*(-(ir**2*yh*(-35*xh**5 + 35*xh**3*(1 + yh**2) - 3*xh*(2 + 5*yh**2))*zh))
        gxy=8.302649259524165E0_r8*(ir**2*(yh**2*(-3 + 5*yh**2) + 5*xh**4*(-1 + 7*yh**2) + xh**2*(3 - 35*yh**4))*zh)
        gxz=8.302649259524165E0_r8*(-(ir**2*yh*(xh**4*(5 - 35*zh**2) + yh**2*(1 - 5*zh**2) + xh**2*(-3 + 15*zh**2 + 5*yh**2*(-1 + 7*zh**2)))))
        gyy=8.302649259524165E0_r8*(ir**2*xh*yh*(-6 + 35*yh**2 - 35*yh**4 + 5*xh**2*(-3 + 7*yh**2))*zh)
        gyz=8.302649259524165E0_r8*(ir**2*xh*(yh**2*(-3 + 15*zh**2 + yh**2*(5 - 35*zh**2)) + xh**2*(1 - 5*zh**2 + 5*yh**2*(-1 + 7*zh**2))))
        gzz=8.302649259524165E0_r8*(5*ir**2*xh*yh*(xh**2 - yh**2)*zh*(-3 + 7*zh**2))
    case(28)
        ! l=5, m=-3, iang=28 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(-1 + 9*zh^2))/16
        ! dYlm/dx = (Sqrt[385/(2*Pi)]*xh*yh*(15*xh^4 - 5*yh^4 + 48*zh^2 + 2*xh^2*(-6 + 5*yh^2 - 60*zh^2) + yh^2*(-4 + 40*zh^2)))/(16*r)
        ! dYlm/dy = (Sqrt[385/(2*Pi)]*(-5*yh^6 + 3*xh^4*(-1 + 5*yh^2) - 24*yh^2*zh^2 + 5*yh^4*(1 + 8*zh^2) + 2*xh^2*(5*yh^4 + 12*zh^2 - 3*yh^2*(1 + 20*zh^2))))/(16*r)
        ! dYlm/dz = -(Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(16 + 5*xh^2 + 5*yh^2 - 40*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (Sqrt[385/(2*Pi)]*yh*(-105*xh^6 - 5*yh^4 + 48*zh^2 + yh^2*(-4 + 40*zh^2) + xh^4*(135 - 70*yh^2 + 840*zh^2) + xh^2*(35*yh^4 + yh^2*(50 - 280*zh^2) - 12*(3 + 50*zh^2))))/(16*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[385/(2*Pi)]*xh*(-35*yh^6 + 15*xh^4*(-1 + 7*yh^2) - 48*zh^2 + 12*yh^2*(1 + 10*zh^2) + 5*yh^4*(1 + 56*zh^2) + 2*xh^2*(6 + 35*yh^4 + 60*zh^2 - 15*yh^2*(3 + 28*zh^2))))/(16*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[385/(2*Pi)]*xh*yh*zh*(105*xh^4 - 35*yh^4 + 10*xh^2*(18 + 7*yh^2 - 84*zh^2) + 48*(-2 + 5*zh^2) + 20*yh^2*(-5 + 14*zh^2)))/(16*r^2)
        ! d^2Ylm/dydx = -(Sqrt[385/(2*Pi)]*xh*(-35*yh^6 + 15*xh^4*(-1 + 7*yh^2) - 48*zh^2 + 12*yh^2*(1 + 10*zh^2) + 5*yh^4*(1 + 56*zh^2) + 2*xh^2*(6 + 35*yh^4 + 60*zh^2 - 15*yh^2*(3 + 28*zh^2))))/(16*r^2)
        ! d^2Ylm/dydy = -(Sqrt[385/(2*Pi)]*yh*(-35*yh^6 + 15*xh^4*(-3 + 7*yh^2) + 48*zh^2 - 20*yh^2*(1 + 14*zh^2) + 5*yh^4*(11 + 56*zh^2) + 2*xh^2*(6 + 35*yh^4 + 180*zh^2 - 35*yh^2*(1 + 12*zh^2))))/(16*r^2)
        ! d^2Ylm/dydz = (Sqrt[385/(2*Pi)]*zh*(35*yh^6 - 15*xh^4*(-1 + 7*yh^2) + yh^4*(55 - 280*zh^2) + 24*yh^2*(-2 + 5*zh^2) - 2*xh^2*(-24 + 35*yh^4 + 60*zh^2 - 105*yh^2*(-1 + 4*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[385/(2*Pi)]*xh*yh*zh*(105*xh^4 - 35*yh^4 + 10*xh^2*(18 + 7*yh^2 - 84*zh^2) + 48*(-2 + 5*zh^2) + 20*yh^2*(-5 + 14*zh^2)))/(16*r^2)
        ! d^2Ylm/dzdy = (Sqrt[385/(2*Pi)]*zh*(35*yh^6 - 15*xh^4*(-1 + 7*yh^2) + yh^4*(55 - 280*zh^2) + 24*yh^2*(-2 + 5*zh^2) - 2*xh^2*(-24 + 35*yh^4 + 60*zh^2 - 105*yh^2*(-1 + 4*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdz = (Sqrt[385/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(5*xh^2*(-1 + 7*zh^2) + 5*yh^2*(-1 + 7*zh^2) - 8*(2 - 25*zh^2 + 35*zh^4)))/(16*r^2)
        y=-4.892382994352504E-1_r8*(yh*(-3*xh**2 + yh**2)*(-1 + 9*zh**2))
        gx=-4.892382994352504E-1_r8*(-(ir*xh*yh*(15*xh**4 - 5*yh**4 + 48*zh**2 + 2*xh**2*(-6 + 5*yh**2 - 60*zh**2) + yh**2*(-4 + 40*zh**2))))
        gy=-4.892382994352504E-1_r8*(-(ir*(-5*yh**6 + 3*xh**4*(-1 + 5*yh**2) - 24*yh**2*zh**2 + 5*yh**4*(1 + 8*zh**2) + 2*xh**2*(5*yh**4 + 12*zh**2 - 3*yh**2*(1 + 20*zh**2)))))
        gz=-4.892382994352504E-1_r8*(ir*yh*(-3*xh**2 + yh**2)*zh*(16 + 5*xh**2 + 5*yh**2 - 40*zh**2))
        gxx=-4.892382994352504E-1_r8*(-(ir**2*yh*(-105*xh**6 - 5*yh**4 + 48*zh**2 + yh**2*(-4 + 40*zh**2) + xh**4*(135 - 70*yh**2 + 840*zh**2) + xh**2*(35*yh**4 + yh**2*(50 - 280*zh**2) - 12*(3 + 50*zh**2)))))
        gxy=-4.892382994352504E-1_r8*(ir**2*xh*(-35*yh**6 + 15*xh**4*(-1 + 7*yh**2) - 48*zh**2 + 12*yh**2*(1 + 10*zh**2) + 5*yh**4*(1 + 56*zh**2) + 2*xh**2*(6 + 35*yh**4 + 60*zh**2 - 15*yh**2*(3 + 28*zh**2))))
        gxz=-4.892382994352504E-1_r8*(ir**2*xh*yh*zh*(105*xh**4 - 35*yh**4 + 10*xh**2*(18 + 7*yh**2 - 84*zh**2) + 48*(-2 + 5*zh**2) + 20*yh**2*(-5 + 14*zh**2)))
        gyy=-4.892382994352504E-1_r8*(ir**2*yh*(-35*yh**6 + 15*xh**4*(-3 + 7*yh**2) + 48*zh**2 - 20*yh**2*(1 + 14*zh**2) + 5*yh**4*(11 + 56*zh**2) + 2*xh**2*(6 + 35*yh**4 + 180*zh**2 - 35*yh**2*(1 + 12*zh**2))))
        gyz=-4.892382994352504E-1_r8*(-(ir**2*zh*(35*yh**6 - 15*xh**4*(-1 + 7*yh**2) + yh**4*(55 - 280*zh**2) + 24*yh**2*(-2 + 5*zh**2) - 2*xh**2*(-24 + 35*yh**4 + 60*zh**2 - 105*yh**2*(-1 + 4*zh**2)))))
        gzz=-4.892382994352504E-1_r8*(-(ir**2*yh*(-3*xh**2 + yh**2)*(5*xh**2*(-1 + 7*zh**2) + 5*yh**2*(-1 + 7*zh**2) - 8*(2 - 25*zh**2 + 35*zh**4))))
    case(29)
        ! l=5, m=-2, iang=29 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*xh*yh*zh*(-1 + 3*zh^2))/4
        ! dYlm/dx = (Sqrt[1155/Pi]*yh*zh*(5*xh^4 - yh^2 + 2*zh^2 + xh^2*(-3 + 5*yh^2 - 10*zh^2)))/(4*r)
        ! dYlm/dy = (Sqrt[1155/Pi]*xh*zh*(5*yh^4 + xh^2*(-1 + 5*yh^2) + 2*zh^2 - yh^2*(3 + 10*zh^2)))/(4*r)
        ! dYlm/dz = (Sqrt[1155/Pi]*xh*yh*(-xh^2 - yh^2 + (6 + 5*xh^2 + 5*yh^2)*zh^2 - 10*zh^4))/(4*r)
        ! d^2Ylm/dxdx = -(Sqrt[1155/Pi]*xh*yh*zh*(6 + 35*xh^4 - 15*yh^2 + 30*zh^2 + 35*xh^2*(-1 + yh^2 - 2*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (Sqrt[1155/Pi]*zh*(5*yh^4 + xh^4*(5 - 35*yh^2) + 2*zh^2 - yh^2*(3 + 10*zh^2) + xh^2*(-3 - 35*yh^4 - 10*zh^2 + 10*yh^2*(3 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[1155/Pi]*yh*(yh^2 - 6*zh^2 - 5*yh^2*zh^2 + 10*zh^4 + 5*xh^4*(-1 + 7*zh^2) + xh^2*(3 + 15*zh^2 - 70*zh^4 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dydx = (Sqrt[1155/Pi]*zh*(5*yh^4 + xh^4*(5 - 35*yh^2) + 2*zh^2 - yh^2*(3 + 10*zh^2) + xh^2*(-3 - 35*yh^4 - 10*zh^2 + 10*yh^2*(3 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dydy = -(Sqrt[1155/Pi]*xh*yh*zh*(6 + 35*yh^4 + 5*xh^2*(-3 + 7*yh^2) + 30*zh^2 - 35*yh^2*(1 + 2*zh^2)))/(4*r^2)
        ! d^2Ylm/dydz = -(Sqrt[1155/Pi]*xh*(2*zh^2*(-3 + 5*zh^2) + 5*yh^4*(-1 + 7*zh^2) + yh^2*(3 + 15*zh^2 - 70*zh^4) + xh^2*(1 - 5*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[1155/Pi]*yh*(yh^2 - 6*zh^2 - 5*yh^2*zh^2 + 10*zh^4 + 5*xh^4*(-1 + 7*zh^2) + xh^2*(3 + 15*zh^2 - 70*zh^4 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdy = -(Sqrt[1155/Pi]*xh*(2*zh^2*(-3 + 5*zh^2) + 5*yh^4*(-1 + 7*zh^2) + yh^2*(3 + 15*zh^2 - 70*zh^4) + xh^2*(1 - 5*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[1155/Pi]*xh*yh*zh*(-12 + 70*zh^2 - 70*zh^4 + 5*xh^2*(-3 + 7*zh^2) + 5*yh^2*(-3 + 7*zh^2)))/(4*r^2)
        y=4.793536784973324E0_r8*(xh*yh*zh*(-1 + 3*zh**2))
        gx=4.793536784973324E0_r8*(ir*yh*zh*(5*xh**4 - yh**2 + 2*zh**2 + xh**2*(-3 + 5*yh**2 - 10*zh**2)))
        gy=4.793536784973324E0_r8*(ir*xh*zh*(5*yh**4 + xh**2*(-1 + 5*yh**2) + 2*zh**2 - yh**2*(3 + 10*zh**2)))
        gz=4.793536784973324E0_r8*(ir*xh*yh*(-xh**2 - yh**2 + (6 + 5*xh**2 + 5*yh**2)*zh**2 - 10*zh**4))
        gxx=4.793536784973324E0_r8*(-(ir**2*xh*yh*zh*(6 + 35*xh**4 - 15*yh**2 + 30*zh**2 + 35*xh**2*(-1 + yh**2 - 2*zh**2))))
        gxy=4.793536784973324E0_r8*(ir**2*zh*(5*yh**4 + xh**4*(5 - 35*yh**2) + 2*zh**2 - yh**2*(3 + 10*zh**2) + xh**2*(-3 - 35*yh**4 - 10*zh**2 + 10*yh**2*(3 + 7*zh**2))))
        gxz=4.793536784973324E0_r8*(-(ir**2*yh*(yh**2 - 6*zh**2 - 5*yh**2*zh**2 + 10*zh**4 + 5*xh**4*(-1 + 7*zh**2) + xh**2*(3 + 15*zh**2 - 70*zh**4 + 5*yh**2*(-1 + 7*zh**2)))))
        gyy=4.793536784973324E0_r8*(-(ir**2*xh*yh*zh*(6 + 35*yh**4 + 5*xh**2*(-3 + 7*yh**2) + 30*zh**2 - 35*yh**2*(1 + 2*zh**2))))
        gyz=4.793536784973324E0_r8*(-(ir**2*xh*(2*zh**2*(-3 + 5*zh**2) + 5*yh**4*(-1 + 7*zh**2) + yh**2*(3 + 15*zh**2 - 70*zh**4) + xh**2*(1 - 5*zh**2 + 5*yh**2*(-1 + 7*zh**2)))))
        gzz=4.793536784973324E0_r8*(-(ir**2*xh*yh*zh*(-12 + 70*zh**2 - 70*zh**4 + 5*xh**2*(-3 + 7*zh**2) + 5*yh**2*(-3 + 7*zh**2))))
    case(30)
        ! l=5, m=-1, iang=30 (=l^2+l+1+m)
        ! Ylm= (Sqrt[165/Pi]*yh*(1 - 14*zh^2 + 21*zh^4))/16
        ! dYlm/dx = -(Sqrt[165/Pi]*xh*yh*(1 - 42*zh^2 + 105*zh^4))/(16*r)
        ! dYlm/dy = -(Sqrt[165/Pi]*(-1 + yh^2 + 14*(1 - 3*yh^2)*zh^2 + 21*(-1 + 5*yh^2)*zh^4))/(16*r)
        ! dYlm/dz = (Sqrt[165/Pi]*yh*zh*(-29 + 21*zh^2*(6 - 5*zh^2)))/(16*r)
        ! d^2Ylm/dxdx = (Sqrt[165/Pi]*yh*(-1 + 42*zh^2 - 105*zh^4 + 3*xh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (Sqrt[165/Pi]*xh*(-1 + 42*zh^2 - 105*zh^4 + 3*yh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[165/Pi]*xh*yh*zh*(29 - 210*zh^2 + 245*zh^4))/(16*r^2)
        ! d^2Ylm/dydx = (Sqrt[165/Pi]*xh*(-1 + 42*zh^2 - 105*zh^4 + 3*yh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[165/Pi]*yh*(-1 + 42*zh^2 - 105*zh^4 + yh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (Sqrt[165/Pi]*zh*(-29 + 126*zh^2 - 105*zh^4 + yh^2*(87 - 630*zh^2 + 735*zh^4)))/(16*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[165/Pi]*xh*yh*zh*(29 - 210*zh^2 + 245*zh^4))/(16*r^2)
        ! d^2Ylm/dzdy = (Sqrt[165/Pi]*zh*(-29 + 126*zh^2 - 105*zh^4 + yh^2*(87 - 630*zh^2 + 735*zh^4)))/(16*r^2)
        ! d^2Ylm/dzdz = (Sqrt[165/Pi]*yh*(-29 + 465*zh^2 - 1155*zh^4 + 735*zh^6))/(16*r^2)
        y=4.529466511956969E-1_r8*(yh*(1 - 14*zh**2 + 21*zh**4))
        gx=4.529466511956969E-1_r8*(-(ir*xh*yh*(1 - 42*zh**2 + 105*zh**4)))
        gy=4.529466511956969E-1_r8*(-(ir*(-1 + yh**2 + 14*(1 - 3*yh**2)*zh**2 + 21*(-1 + 5*yh**2)*zh**4)))
        gz=4.529466511956969E-1_r8*(ir*yh*zh*(-29 + 21*zh**2*(6 - 5*zh**2)))
        gxx=4.529466511956969E-1_r8*(ir**2*yh*(-1 + 42*zh**2 - 105*zh**4 + 3*xh**2*(1 - 70*zh**2 + 245*zh**4)))
        gxy=4.529466511956969E-1_r8*(ir**2*xh*(-1 + 42*zh**2 - 105*zh**4 + 3*yh**2*(1 - 70*zh**2 + 245*zh**4)))
        gxz=4.529466511956969E-1_r8*(3*ir**2*xh*yh*zh*(29 - 210*zh**2 + 245*zh**4))
        gyy=4.529466511956969E-1_r8*(3*ir**2*yh*(-1 + 42*zh**2 - 105*zh**4 + yh**2*(1 - 70*zh**2 + 245*zh**4)))
        gyz=4.529466511956969E-1_r8*(ir**2*zh*(-29 + 126*zh**2 - 105*zh**4 + yh**2*(87 - 630*zh**2 + 735*zh**4)))
        gzz=4.529466511956969E-1_r8*(ir**2*yh*(-29 + 465*zh**2 - 1155*zh**4 + 735*zh**6))
    case(31)
        ! l=5, m=0, iang=31 (=l^2+l+1+m)
        ! Ylm= (Sqrt[11/Pi]*zh*(15 - 70*zh^2 + 63*zh^4))/16
        ! dYlm/dx = (-15*Sqrt[11/Pi]*xh*zh*(1 - 14*zh^2 + 21*zh^4))/(16*r)
        ! dYlm/dy = (-15*Sqrt[11/Pi]*yh*zh*(1 - 14*zh^2 + 21*zh^4))/(16*r)
        ! dYlm/dz = (15*Sqrt[11/Pi]*(1 - 15*zh^2 + 35*zh^4 - 21*zh^6))/(16*r)
        ! d^2Ylm/dxdx = (15*Sqrt[11/Pi]*zh*(-1 + 14*zh^2 - 21*zh^4 + xh^2*(3 - 70*zh^2 + 147*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (15*Sqrt[11/Pi]*xh*yh*zh*(3 - 70*zh^2 + 147*zh^4))/(16*r^2)
        ! d^2Ylm/dxdz = (15*Sqrt[11/Pi]*xh*(-1 + 45*zh^2 - 175*zh^4 + 147*zh^6))/(16*r^2)
        ! d^2Ylm/dydx = (15*Sqrt[11/Pi]*xh*yh*zh*(3 - 70*zh^2 + 147*zh^4))/(16*r^2)
        ! d^2Ylm/dydy = (15*Sqrt[11/Pi]*zh*(-1 + 14*zh^2 - 21*zh^4 + yh^2*(3 - 70*zh^2 + 147*zh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (15*Sqrt[11/Pi]*yh*(-1 + 45*zh^2 - 175*zh^4 + 147*zh^6))/(16*r^2)
        ! d^2Ylm/dzdx = (15*Sqrt[11/Pi]*xh*(-1 + 45*zh^2 - 175*zh^4 + 147*zh^6))/(16*r^2)
        ! d^2Ylm/dzdy = (15*Sqrt[11/Pi]*yh*(-1 + 45*zh^2 - 175*zh^4 + 147*zh^6))/(16*r^2)
        ! d^2Ylm/dzdz = (15*Sqrt[11/Pi]*zh*(-31 + 185*zh^2 - 301*zh^4 + 147*zh^6))/(16*r^2)
        y=1.169503224534236E-1_r8*(zh*(15 - 70*zh**2 + 63*zh**4))
        gx=1.169503224534236E-1_r8*(-15*ir*xh*zh*(1 - 14*zh**2 + 21*zh**4))
        gy=1.169503224534236E-1_r8*(-15*ir*yh*zh*(1 - 14*zh**2 + 21*zh**4))
        gz=1.169503224534236E-1_r8*(15*ir*(1 - 15*zh**2 + 35*zh**4 - 21*zh**6))
        gxx=1.169503224534236E-1_r8*(15*ir**2*zh*(-1 + 14*zh**2 - 21*zh**4 + xh**2*(3 - 70*zh**2 + 147*zh**4)))
        gxy=1.169503224534236E-1_r8*(15*ir**2*xh*yh*zh*(3 - 70*zh**2 + 147*zh**4))
        gxz=1.169503224534236E-1_r8*(15*ir**2*xh*(-1 + 45*zh**2 - 175*zh**4 + 147*zh**6))
        gyy=1.169503224534236E-1_r8*(15*ir**2*zh*(-1 + 14*zh**2 - 21*zh**4 + yh**2*(3 - 70*zh**2 + 147*zh**4)))
        gyz=1.169503224534236E-1_r8*(15*ir**2*yh*(-1 + 45*zh**2 - 175*zh**4 + 147*zh**6))
        gzz=1.169503224534236E-1_r8*(15*ir**2*zh*(-31 + 185*zh**2 - 301*zh**4 + 147*zh**6))
    case(32)
        ! l=5, m=1, iang=32 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[165/Pi]*xh*(1 - 14*zh^2 + 21*zh^4))/16
        ! dYlm/dx = (Sqrt[165/Pi]*(-1 + xh^2 + 14*(1 - 3*xh^2)*zh^2 + 21*(-1 + 5*xh^2)*zh^4))/(16*r)
        ! dYlm/dy = (Sqrt[165/Pi]*xh*yh*(1 - 42*zh^2 + 105*zh^4))/(16*r)
        ! dYlm/dz = (Sqrt[165/Pi]*xh*zh*(29 + 21*zh^2*(-6 + 5*zh^2)))/(16*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[165/Pi]*xh*(-1 + 42*zh^2 - 105*zh^4 + xh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[165/Pi]*yh*(-1 + 42*zh^2 - 105*zh^4 + 3*xh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[165/Pi]*zh*(-29 + 126*zh^2 - 105*zh^4 + xh^2*(87 - 630*zh^2 + 735*zh^4)))/(16*r^2)
        ! d^2Ylm/dydx = -(Sqrt[165/Pi]*yh*(-1 + 42*zh^2 - 105*zh^4 + 3*xh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dydy = -(Sqrt[165/Pi]*xh*(-1 + 42*zh^2 - 105*zh^4 + 3*yh^2*(1 - 70*zh^2 + 245*zh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[165/Pi]*xh*yh*zh*(29 - 210*zh^2 + 245*zh^4))/(16*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[165/Pi]*zh*(-29 + 126*zh^2 - 105*zh^4 + xh^2*(87 - 630*zh^2 + 735*zh^4)))/(16*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[165/Pi]*xh*yh*zh*(29 - 210*zh^2 + 245*zh^4))/(16*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[165/Pi]*xh*(-29 + 465*zh^2 - 1155*zh^4 + 735*zh^6))/(16*r^2)
        y=-4.529466511956969E-1_r8*(xh*(1 - 14*zh**2 + 21*zh**4))
        gx=-4.529466511956969E-1_r8*(-(ir*(-1 + xh**2 + 14*(1 - 3*xh**2)*zh**2 + 21*(-1 + 5*xh**2)*zh**4)))
        gy=-4.529466511956969E-1_r8*(-(ir*xh*yh*(1 - 42*zh**2 + 105*zh**4)))
        gz=-4.529466511956969E-1_r8*(-(ir*xh*zh*(29 + 21*zh**2*(-6 + 5*zh**2))))
        gxx=-4.529466511956969E-1_r8*(3*ir**2*xh*(-1 + 42*zh**2 - 105*zh**4 + xh**2*(1 - 70*zh**2 + 245*zh**4)))
        gxy=-4.529466511956969E-1_r8*(ir**2*yh*(-1 + 42*zh**2 - 105*zh**4 + 3*xh**2*(1 - 70*zh**2 + 245*zh**4)))
        gxz=-4.529466511956969E-1_r8*(ir**2*zh*(-29 + 126*zh**2 - 105*zh**4 + xh**2*(87 - 630*zh**2 + 735*zh**4)))
        gyy=-4.529466511956969E-1_r8*(ir**2*xh*(-1 + 42*zh**2 - 105*zh**4 + 3*yh**2*(1 - 70*zh**2 + 245*zh**4)))
        gyz=-4.529466511956969E-1_r8*(3*ir**2*xh*yh*zh*(29 - 210*zh**2 + 245*zh**4))
        gzz=-4.529466511956969E-1_r8*(ir**2*xh*(-29 + 465*zh**2 - 1155*zh**4 + 735*zh**6))
    case(33)
        ! l=5, m=2, iang=33 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1155/Pi]*(xh^2 - yh^2)*zh*(-1 + 3*zh^2))/8
        ! dYlm/dx = (Sqrt[1155/Pi]*xh*zh*(-4*xh^2 + 5*xh^4 - 5*yh^4 + 2*(2 - 5*xh^2 + 5*yh^2)*zh^2))/(8*r)
        ! dYlm/dy = (Sqrt[1155/Pi]*yh*zh*(5*xh^4 + 4*yh^2 - 5*yh^4 - 2*(2 + 5*xh^2 - 5*yh^2)*zh^2))/(8*r)
        ! dYlm/dz = (Sqrt[1155/Pi]*(xh - yh)*(xh + yh)*(-xh^2 - yh^2 + (6 + 5*xh^2 + 5*yh^2)*zh^2 - 10*zh^4))/(8*r)
        ! d^2Ylm/dxdx = (Sqrt[1155/Pi]*zh*(-35*xh^6 - 5*yh^4 + 4*zh^2 + 10*yh^2*zh^2 + 5*xh^4*(9 + 14*zh^2) + xh^2*(-12 + 35*yh^4 - 50*zh^2 - 70*yh^2*zh^2)))/(8*r^2)
        ! d^2Ylm/dxdy = (-5*Sqrt[1155/Pi]*xh*yh*(xh^2 - yh^2)*zh*(-4 + 7*xh^2 + 7*yh^2 - 14*zh^2))/(8*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[1155/Pi]*xh*(yh^4*(5 - 35*zh^2) + 4*zh^2*(-3 + 5*zh^2) + 10*yh^2*zh^2*(-3 + 7*zh^2) + 5*xh^4*(-1 + 7*zh^2) + xh^2*(4 + 10*zh^2 - 70*zh^4)))/(8*r^2)
        ! d^2Ylm/dydx = (-5*Sqrt[1155/Pi]*xh*yh*(xh^2 - yh^2)*zh*(-4 + 7*xh^2 + 7*yh^2 - 14*zh^2))/(8*r^2)
        ! d^2Ylm/dydy = (Sqrt[1155/Pi]*zh*(35*yh^6 + xh^4*(5 - 35*yh^2) - 4*zh^2 + 10*xh^2*(-1 + 7*yh^2)*zh^2 - 5*yh^4*(9 + 14*zh^2) + 2*yh^2*(6 + 25*zh^2)))/(8*r^2)
        ! d^2Ylm/dydz = (Sqrt[1155/Pi]*yh*(xh^4*(5 - 35*zh^2) + 4*zh^2*(-3 + 5*zh^2) + 10*xh^2*zh^2*(-3 + 7*zh^2) + 5*yh^4*(-1 + 7*zh^2) + yh^2*(4 + 10*zh^2 - 70*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[1155/Pi]*xh*(yh^4*(5 - 35*zh^2) + 4*zh^2*(-3 + 5*zh^2) + 10*yh^2*zh^2*(-3 + 7*zh^2) + 5*xh^4*(-1 + 7*zh^2) + xh^2*(4 + 10*zh^2 - 70*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdy = (Sqrt[1155/Pi]*yh*(xh^4*(5 - 35*zh^2) + 4*zh^2*(-3 + 5*zh^2) + 10*xh^2*zh^2*(-3 + 7*zh^2) + 5*yh^4*(-1 + 7*zh^2) + yh^2*(4 + 10*zh^2 - 70*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[1155/Pi]*(xh^2 - yh^2)*zh*(-12 + 70*zh^2 - 70*zh^4 + 5*xh^2*(-3 + 7*zh^2) + 5*yh^2*(-3 + 7*zh^2)))/(8*r^2)
        y=2.396768392486662E0_r8*((xh**2 - yh**2)*zh*(-1 + 3*zh**2))
        gx=2.396768392486662E0_r8*(ir*xh*zh*(-4*xh**2 + 5*xh**4 - 5*yh**4 + 2*(2 - 5*xh**2 + 5*yh**2)*zh**2))
        gy=2.396768392486662E0_r8*(ir*yh*zh*(5*xh**4 + 4*yh**2 - 5*yh**4 - 2*(2 + 5*xh**2 - 5*yh**2)*zh**2))
        gz=2.396768392486662E0_r8*(ir*(xh - yh)*(xh + yh)*(-xh**2 - yh**2 + (6 + 5*xh**2 + 5*yh**2)*zh**2 - 10*zh**4))
        gxx=2.396768392486662E0_r8*(ir**2*zh*(-35*xh**6 - 5*yh**4 + 4*zh**2 + 10*yh**2*zh**2 + 5*xh**4*(9 + 14*zh**2) + xh**2*(-12 + 35*yh**4 - 50*zh**2 - 70*yh**2*zh**2)))
        gxy=2.396768392486662E0_r8*(-5*ir**2*xh*yh*(xh**2 - yh**2)*zh*(-4 + 7*xh**2 + 7*yh**2 - 14*zh**2))
        gxz=2.396768392486662E0_r8*(-(ir**2*xh*(yh**4*(5 - 35*zh**2) + 4*zh**2*(-3 + 5*zh**2) + 10*yh**2*zh**2*(-3 + 7*zh**2) + 5*xh**4*(-1 + 7*zh**2) + xh**2*(4 + 10*zh**2 - 70*zh**4))))
        gyy=2.396768392486662E0_r8*(ir**2*zh*(35*yh**6 + xh**4*(5 - 35*yh**2) - 4*zh**2 + 10*xh**2*(-1 + 7*yh**2)*zh**2 - 5*yh**4*(9 + 14*zh**2) + 2*yh**2*(6 + 25*zh**2)))
        gyz=2.396768392486662E0_r8*(ir**2*yh*(xh**4*(5 - 35*zh**2) + 4*zh**2*(-3 + 5*zh**2) + 10*xh**2*zh**2*(-3 + 7*zh**2) + 5*yh**4*(-1 + 7*zh**2) + yh**2*(4 + 10*zh**2 - 70*zh**4)))
        gzz=2.396768392486662E0_r8*(-(ir**2*(xh**2 - yh**2)*zh*(-12 + 70*zh**2 - 70*zh**4 + 5*xh**2*(-3 + 7*zh**2) + 5*yh**2*(-3 + 7*zh**2))))
    case(34)
        ! l=5, m=3, iang=34 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(-1 + 9*zh^2))/16
        ! dYlm/dx = -(Sqrt[385/(2*Pi)]*(5*xh^6 - 5*xh^4*(1 + 2*yh^2 + 8*zh^2) + 3*(yh^4 - 8*yh^2*zh^2) + 3*xh^2*(-5*yh^4 + 8*zh^2 + yh^2*(2 + 40*zh^2))))/(16*r)
        ! dYlm/dy = (Sqrt[385/(2*Pi)]*xh*yh*(-5*xh^4 + 48*zh^2 + 3*yh^2*(-4 + 5*yh^2 - 40*zh^2) + 2*xh^2*(-2 + 5*yh^2 + 20*zh^2)))/(16*r)
        ! dYlm/dz = -(Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(16 + 5*xh^2 + 5*yh^2 - 40*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (Sqrt[385/(2*Pi)]*xh*(35*xh^6 + 45*yh^4 - 48*zh^2 - 12*yh^2*(1 + 30*zh^2) - 5*xh^4*(11 + 14*yh^2 + 56*zh^2) + 5*xh^2*(4 - 21*yh^4 + 56*zh^2 + 14*yh^2*(1 + 12*zh^2))))/(16*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[385/(2*Pi)]*yh*(-35*xh^6 + 5*xh^4*(1 + 14*yh^2 + 56*zh^2) + 3*xh^2*(4 + 35*yh^4 + 40*zh^2 - 10*yh^2*(3 + 28*zh^2)) + 3*(-5*yh^4 - 16*zh^2 + yh^2*(4 + 40*zh^2))))/(16*r^2)
        ! d^2Ylm/dxdz = -(Sqrt[385/(2*Pi)]*zh*(-35*xh^6 - 3*yh^2*(16 + 5*yh^2 - 40*zh^2) + 5*xh^4*(-11 + 14*yh^2 + 56*zh^2) + 3*xh^2*(16 + 35*yh^4 - 40*zh^2 - 70*yh^2*(-1 + 4*zh^2))))/(16*r^2)
        ! d^2Ylm/dydx = -(Sqrt[385/(2*Pi)]*yh*(-35*xh^6 + 5*xh^4*(1 + 14*yh^2 + 56*zh^2) + 3*xh^2*(4 + 35*yh^4 + 40*zh^2 - 10*yh^2*(3 + 28*zh^2)) + 3*(-5*yh^4 - 16*zh^2 + yh^2*(4 + 40*zh^2))))/(16*r^2)
        ! d^2Ylm/dydy = (Sqrt[385/(2*Pi)]*xh*(5*xh^4*(-1 + 7*yh^2) + xh^2*(-4 - 70*yh^4 + 40*zh^2 + yh^2*(50 - 280*zh^2)) + 3*(-35*yh^6 + 16*zh^2 - 4*yh^2*(3 + 50*zh^2) + 5*yh^4*(9 + 56*zh^2))))/(16*r^2)
        ! d^2Ylm/dydz = (Sqrt[385/(2*Pi)]*xh*yh*zh*(35*xh^4 - 10*xh^2*(-10 + 7*yh^2 + 28*zh^2) - 3*(-32 + 35*yh^4 + 80*zh^2 + yh^2*(60 - 280*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdx = -(Sqrt[385/(2*Pi)]*zh*(-35*xh^6 - 3*yh^2*(16 + 5*yh^2 - 40*zh^2) + 5*xh^4*(-11 + 14*yh^2 + 56*zh^2) + 3*xh^2*(16 + 35*yh^4 - 40*zh^2 - 70*yh^2*(-1 + 4*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdy = (Sqrt[385/(2*Pi)]*xh*yh*zh*(35*xh^4 - 10*xh^2*(-10 + 7*yh^2 + 28*zh^2) - 3*(-32 + 35*yh^4 + 80*zh^2 + yh^2*(60 - 280*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdz = (Sqrt[385/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(5*xh^2*(-1 + 7*zh^2) + 5*yh^2*(-1 + 7*zh^2) - 8*(2 - 25*zh^2 + 35*zh^4)))/(16*r^2)
        y=-4.892382994352504E-1_r8*(xh*(xh**2 - 3*yh**2)*(-1 + 9*zh**2))
        gx=-4.892382994352504E-1_r8*(ir*(5*xh**6 - 5*xh**4*(1 + 2*yh**2 + 8*zh**2) + 3*(yh**4 - 8*yh**2*zh**2) + 3*xh**2*(-5*yh**4 + 8*zh**2 + yh**2*(2 + 40*zh**2))))
        gy=-4.892382994352504E-1_r8*(-(ir*xh*yh*(-5*xh**4 + 48*zh**2 + 3*yh**2*(-4 + 5*yh**2 - 40*zh**2) + 2*xh**2*(-2 + 5*yh**2 + 20*zh**2))))
        gz=-4.892382994352504E-1_r8*(ir*xh*(xh**2 - 3*yh**2)*zh*(16 + 5*xh**2 + 5*yh**2 - 40*zh**2))
        gxx=-4.892382994352504E-1_r8*(-(ir**2*xh*(35*xh**6 + 45*yh**4 - 48*zh**2 - 12*yh**2*(1 + 30*zh**2) - 5*xh**4*(11 + 14*yh**2 + 56*zh**2) + 5*xh**2*(4 - 21*yh**4 + 56*zh**2 + 14*yh**2*(1 + 12*zh**2)))))
        gxy=-4.892382994352504E-1_r8*(ir**2*yh*(-35*xh**6 + 5*xh**4*(1 + 14*yh**2 + 56*zh**2) + 3*xh**2*(4 + 35*yh**4 + 40*zh**2 - 10*yh**2*(3 + 28*zh**2)) + 3*(-5*yh**4 - 16*zh**2 + yh**2*(4 + 40*zh**2))))
        gxz=-4.892382994352504E-1_r8*(ir**2*zh*(-35*xh**6 - 3*yh**2*(16 + 5*yh**2 - 40*zh**2) + 5*xh**4*(-11 + 14*yh**2 + 56*zh**2) + 3*xh**2*(16 + 35*yh**4 - 40*zh**2 - 70*yh**2*(-1 + 4*zh**2))))
        gyy=-4.892382994352504E-1_r8*(-(ir**2*xh*(5*xh**4*(-1 + 7*yh**2) + xh**2*(-4 - 70*yh**4 + 40*zh**2 + yh**2*(50 - 280*zh**2)) + 3*(-35*yh**6 + 16*zh**2 - 4*yh**2*(3 + 50*zh**2) + 5*yh**4*(9 + 56*zh**2)))))
        gyz=-4.892382994352504E-1_r8*(-(ir**2*xh*yh*zh*(35*xh**4 - 10*xh**2*(-10 + 7*yh**2 + 28*zh**2) - 3*(-32 + 35*yh**4 + 80*zh**2 + yh**2*(60 - 280*zh**2)))))
        gzz=-4.892382994352504E-1_r8*(-(ir**2*xh*(xh**2 - 3*yh**2)*(5*xh**2*(-1 + 7*zh**2) + 5*yh**2*(-1 + 7*zh**2) - 8*(2 - 25*zh**2 + 35*zh**4))))
    case(35)
        ! l=5, m=4, iang=35 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh)/16
        ! dYlm/dx = (-3*Sqrt[385/Pi]*xh*(xh^2*(-4 + 5*xh^2) + 6*(2 - 5*xh^2)*yh^2 + 5*yh^4)*zh)/(16*r)
        ! dYlm/dy = (-3*Sqrt[385/Pi]*yh*(5*xh^4 - 4*yh^2 + 5*yh^4 + 6*xh^2*(2 - 5*yh^2))*zh)/(16*r)
        ! dYlm/dz = (-3*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 5*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (3*Sqrt[385/Pi]*(35*xh^6 - yh^2*(12 + 5*yh^2) - 15*xh^4*(3 + 14*yh^2) + xh^2*(12 + 150*yh^2 + 35*yh^4))*zh)/(16*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[385/Pi]*xh*yh*(-24 + 35*xh^4 + 40*yh^2 + 35*yh^4 + xh^2*(40 - 210*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[385/Pi]*(5*xh^5*(-1 + 7*zh^2) + xh^3*(4 - 20*zh^2 - 30*yh^2*(-1 + 7*zh^2)) + xh*yh^2*(-12 + 60*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(16*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[385/Pi]*xh*yh*(-24 + 35*xh^4 + 40*yh^2 + 35*yh^4 + xh^2*(40 - 210*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[385/Pi]*(5*xh^4*(-1 + 7*yh^2) + yh^2*(12 - 45*yh^2 + 35*yh^4) - 6*xh^2*(2 - 25*yh^2 + 35*yh^4))*zh)/(16*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[385/Pi]*yh*(5*xh^4*(-1 + 7*zh^2) + yh^2*(4 - 20*zh^2 + 5*yh^2*(-1 + 7*zh^2)) - 6*xh^2*(2 - 10*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[385/Pi]*(5*xh^5*(-1 + 7*zh^2) + xh^3*(4 - 20*zh^2 - 30*yh^2*(-1 + 7*zh^2)) + xh*yh^2*(-12 + 60*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[385/Pi]*yh*(5*xh^4*(-1 + 7*zh^2) + yh^2*(4 - 20*zh^2 + 5*yh^2*(-1 + 7*zh^2)) - 6*xh^2*(2 - 10*zh^2 + 5*yh^2*(-1 + 7*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdz = (15*Sqrt[385/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh*(-3 + 7*zh^2))/(16*r^2)
        y=2.075662314881041E0_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*zh)
        gx=2.075662314881041E0_r8*(-(ir*xh*(xh**2*(-4 + 5*xh**2) + 6*(2 - 5*xh**2)*yh**2 + 5*yh**4)*zh))
        gy=2.075662314881041E0_r8*(-(ir*yh*(5*xh**4 - 4*yh**2 + 5*yh**4 + 6*xh**2*(2 - 5*yh**2))*zh))
        gz=2.075662314881041E0_r8*(-(ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 5*zh**2)))
        gxx=2.075662314881041E0_r8*(ir**2*(35*xh**6 - yh**2*(12 + 5*yh**2) - 15*xh**4*(3 + 14*yh**2) + xh**2*(12 + 150*yh**2 + 35*yh**4))*zh)
        gxy=2.075662314881041E0_r8*(ir**2*xh*yh*(-24 + 35*xh**4 + 40*yh**2 + 35*yh**4 + xh**2*(40 - 210*yh**2))*zh)
        gxz=2.075662314881041E0_r8*(ir**2*(5*xh**5*(-1 + 7*zh**2) + xh**3*(4 - 20*zh**2 - 30*yh**2*(-1 + 7*zh**2)) + xh*yh**2*(-12 + 60*zh**2 + 5*yh**2*(-1 + 7*zh**2))))
        gyy=2.075662314881041E0_r8*(ir**2*(5*xh**4*(-1 + 7*yh**2) + yh**2*(12 - 45*yh**2 + 35*yh**4) - 6*xh**2*(2 - 25*yh**2 + 35*yh**4))*zh)
        gyz=2.075662314881041E0_r8*(ir**2*yh*(5*xh**4*(-1 + 7*zh**2) + yh**2*(4 - 20*zh**2 + 5*yh**2*(-1 + 7*zh**2)) - 6*xh**2*(2 - 10*zh**2 + 5*yh**2*(-1 + 7*zh**2))))
        gzz=2.075662314881041E0_r8*(5*ir**2*(xh**4 - 6*xh**2*yh**2 + yh**4)*zh*(-3 + 7*zh**2))
    case(36)
        ! l=5, m=5, iang=36 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4))/16
        ! dYlm/dx = (15*Sqrt[77/(2*Pi)]*(xh^6 - yh^4 + xh^2*yh^2*(6 + 5*yh^2) - xh^4*(1 + 10*yh^2)))/(16*r)
        ! dYlm/dy = (15*Sqrt[77/(2*Pi)]*xh*yh*(xh^2*(4 + xh^2) - 2*(2 + 5*xh^2)*yh^2 + 5*yh^4))/(16*r)
        ! dYlm/dz = (15*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh)/(16*r)
        ! d^2Ylm/dxdx = (-15*Sqrt[77/(2*Pi)]*(7*xh^7 - 3*xh*yh^2*(4 + 5*yh^2) - xh^5*(11 + 70*yh^2) + xh^3*(4 + 70*yh^2 + 35*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (-15*Sqrt[77/(2*Pi)]*yh*(7*xh^6 + 4*yh^2 - 5*yh^4 + xh^4*(15 - 70*yh^2) + xh^2*(-12 + 10*yh^2 + 35*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdz = (-15*Sqrt[77/(2*Pi)]*(7*xh^6 - 5*yh^4 + 5*xh^2*yh^2*(6 + 7*yh^2) - 5*xh^4*(1 + 14*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dydx = (-15*Sqrt[77/(2*Pi)]*yh*(7*xh^6 + 4*yh^2 - 5*yh^4 + xh^4*(15 - 70*yh^2) + xh^2*(-12 + 10*yh^2 + 35*yh^4)))/(16*r^2)
        ! d^2Ylm/dydy = (-15*Sqrt[77/(2*Pi)]*xh*(xh^4*(-1 + 7*yh^2) + xh^2*(-4 + 50*yh^2 - 70*yh^4) + yh^2*(12 - 45*yh^2 + 35*yh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (-15*Sqrt[77/(2*Pi)]*xh*yh*(7*xh^4 + xh^2*(20 - 70*yh^2) + 5*yh^2*(-4 + 7*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dzdx = (-15*Sqrt[77/(2*Pi)]*(7*xh^6 - 5*yh^4 + 5*xh^2*yh^2*(6 + 7*yh^2) - 5*xh^4*(1 + 14*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dzdy = (-15*Sqrt[77/(2*Pi)]*xh*yh*(7*xh^4 + xh^2*(20 - 70*yh^2) + 5*yh^2*(-4 + 7*yh^2))*zh)/(16*r^2)
        ! d^2Ylm/dzdz = (-15*Sqrt[77/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*(-1 + 7*zh^2))/(16*r^2)
        y=-6.563820568401701E-1_r8*(xh**5 - 10*xh**3*yh**2 + 5*xh*yh**4)
        gx=-6.563820568401701E-1_r8*(-5*ir*(xh**6 - yh**4 + xh**2*yh**2*(6 + 5*yh**2) - xh**4*(1 + 10*yh**2)))
        gy=-6.563820568401701E-1_r8*(-5*ir*xh*yh*(xh**2*(4 + xh**2) - 2*(2 + 5*xh**2)*yh**2 + 5*yh**4))
        gz=-6.563820568401701E-1_r8*(-5*ir*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh)
        gxx=-6.563820568401701E-1_r8*(5*ir**2*(7*xh**7 - 3*xh*yh**2*(4 + 5*yh**2) - xh**5*(11 + 70*yh**2) + xh**3*(4 + 70*yh**2 + 35*yh**4)))
        gxy=-6.563820568401701E-1_r8*(5*ir**2*yh*(7*xh**6 + 4*yh**2 - 5*yh**4 + xh**4*(15 - 70*yh**2) + xh**2*(-12 + 10*yh**2 + 35*yh**4)))
        gxz=-6.563820568401701E-1_r8*(5*ir**2*(7*xh**6 - 5*yh**4 + 5*xh**2*yh**2*(6 + 7*yh**2) - 5*xh**4*(1 + 14*yh**2))*zh)
        gyy=-6.563820568401701E-1_r8*(5*ir**2*xh*(xh**4*(-1 + 7*yh**2) + xh**2*(-4 + 50*yh**2 - 70*yh**4) + yh**2*(12 - 45*yh**2 + 35*yh**4)))
        gyz=-6.563820568401701E-1_r8*(5*ir**2*xh*yh*(7*xh**4 + xh**2*(20 - 70*yh**2) + 5*yh**2*(-4 + 7*yh**2))*zh)
        gzz=-6.563820568401701E-1_r8*(5*ir**2*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*(-1 + 7*zh**2))
    case(37)
        ! l=6, m=-6, iang=37 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4))/16
        ! dYlm/dx = (3*Sqrt[3003/(2*Pi)]*yh*(-6*xh^6 + yh^4 - 2*xh^2*yh^2*(5 + 3*yh^2) + 5*xh^4*(1 + 4*yh^2)))/(16*r)
        ! dYlm/dy = (3*Sqrt[3003/(2*Pi)]*xh*(xh^4 - 2*xh^2*(5 + 3*xh^2)*yh^2 + 5*(1 + 4*xh^2)*yh^4 - 6*yh^6))/(16*r)
        ! dYlm/dz = (-3*Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4)*zh)/(8*r)
        ! d^2Ylm/dxdx = (3*Sqrt[3003/(2*Pi)]*yh*(24*xh^7 - xh*yh^2*(10 + 9*yh^2) - xh^5*(33 + 80*yh^2) + 2*xh^3*(5 + 35*yh^2 + 12*yh^4)))/(8*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[3003/(2*Pi)]*(5*yh^4 - 6*yh^6 + 6*xh^6*(-1 + 8*yh^2) + xh^4*(5 + 30*yh^2 - 160*yh^4) + 6*xh^2*yh^2*(-5 + 5*yh^2 + 8*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[3003/(2*Pi)]*yh*(24*xh^6 - 3*yh^4 + 6*xh^2*yh^2*(5 + 4*yh^2) - 5*xh^4*(3 + 16*yh^2))*zh)/(8*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[3003/(2*Pi)]*(5*yh^4 - 6*yh^6 + 6*xh^6*(-1 + 8*yh^2) + xh^4*(5 + 30*yh^2 - 160*yh^4) + 6*xh^2*yh^2*(-5 + 5*yh^2 + 8*yh^4)))/(16*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4*(-3 + 8*yh^2) - 10*xh^2*(1 - 7*yh^2 + 8*yh^4) + yh^2*(10 - 33*yh^2 + 24*yh^4)))/(8*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[3003/(2*Pi)]*xh*(3*yh^4*(-5 + 8*yh^2) + 3*xh^4*(-1 + 8*yh^2) + xh^2*(30*yh^2 - 80*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[3003/(2*Pi)]*yh*(24*xh^6 - 3*yh^4 + 6*xh^2*yh^2*(5 + 4*yh^2) - 5*xh^4*(3 + 16*yh^2))*zh)/(8*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[3003/(2*Pi)]*xh*(3*yh^4*(-5 + 8*yh^2) + 3*xh^4*(-1 + 8*yh^2) + xh^2*(30*yh^2 - 80*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[3003/(2*Pi)]*xh*yh*(3*xh^4 - 10*xh^2*yh^2 + 3*yh^4)*(-1 + 8*zh^2))/(8*r^2)
        y=1.366368210383829E0_r8*(3*xh**5*yh - 10*xh**3*yh**3 + 3*xh*yh**5)
        gx=1.366368210383829E0_r8*(3*ir*yh*(-6*xh**6 + yh**4 - 2*xh**2*yh**2*(5 + 3*yh**2) + 5*xh**4*(1 + 4*yh**2)))
        gy=1.366368210383829E0_r8*(3*ir*xh*(xh**4 - 2*xh**2*(5 + 3*xh**2)*yh**2 + 5*(1 + 4*xh**2)*yh**4 - 6*yh**6))
        gz=1.366368210383829E0_r8*(-6*ir*xh*yh*(3*xh**4 - 10*xh**2*yh**2 + 3*yh**4)*zh)
        gxx=1.366368210383829E0_r8*(6*ir**2*yh*(24*xh**7 - xh*yh**2*(10 + 9*yh**2) - xh**5*(33 + 80*yh**2) + 2*xh**3*(5 + 35*yh**2 + 12*yh**4)))
        gxy=1.366368210383829E0_r8*(3*ir**2*(5*yh**4 - 6*yh**6 + 6*xh**6*(-1 + 8*yh**2) + xh**4*(5 + 30*yh**2 - 160*yh**4) + 6*xh**2*yh**2*(-5 + 5*yh**2 + 8*yh**4)))
        gxz=1.366368210383829E0_r8*(6*ir**2*yh*(24*xh**6 - 3*yh**4 + 6*xh**2*yh**2*(5 + 4*yh**2) - 5*xh**4*(3 + 16*yh**2))*zh)
        gyy=1.366368210383829E0_r8*(6*ir**2*xh*yh*(3*xh**4*(-3 + 8*yh**2) - 10*xh**2*(1 - 7*yh**2 + 8*yh**4) + yh**2*(10 - 33*yh**2 + 24*yh**4)))
        gyz=1.366368210383829E0_r8*(6*ir**2*xh*(3*yh**4*(-5 + 8*yh**2) + 3*xh**4*(-1 + 8*yh**2) + xh**2*(30*yh**2 - 80*yh**4))*zh)
        gzz=1.366368210383829E0_r8*(6*ir**2*xh*yh*(3*xh**4 - 10*xh**2*yh**2 + 3*yh**4)*(-1 + 8*zh**2))
    case(38)
        ! l=6, m=-5, iang=38 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh)/16
        ! dYlm/dx = (-3*Sqrt[1001/(2*Pi)]*xh*yh*(15*xh^4 + 10*yh^2 + 3*yh^4 - 10*xh^2*(1 + 3*yh^2))*zh)/(8*r)
        ! dYlm/dy = (-3*Sqrt[1001/(2*Pi)]*(-5*xh^4 + 30*xh^2*(1 + xh^2)*yh^2 - 5*(1 + 12*xh^2)*yh^4 + 6*yh^6)*zh)/(16*r)
        ! dYlm/dz = (-3*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*(-1 + 6*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (3*Sqrt[1001/(2*Pi)]*yh*(120*xh^6 - yh^2*(10 + 3*yh^2) - 15*xh^4*(9 + 16*yh^2) + 6*xh^2*(5 + 25*yh^2 + 4*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[1001/(2*Pi)]*(15*xh^5*(-1 + 8*yh^2) + xh^3*(10 + 30*yh^2 - 240*yh^4) + 3*xh*yh^2*(-10 + 15*yh^2 + 8*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[1001/(2*Pi)]*yh*(15*xh^5*(-1 + 8*zh^2) - 10*xh^3*(-1 + 6*zh^2 + 3*yh^2*(-1 + 8*zh^2)) + xh*yh^2*(-10 + 60*zh^2 + 3*yh^2*(-1 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[1001/(2*Pi)]*(15*xh^5*(-1 + 8*yh^2) + xh^3*(10 + 30*yh^2 - 240*yh^4) + 3*xh*yh^2*(-10 + 15*yh^2 + 8*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[1001/(2*Pi)]*yh*(15*xh^4*(-3 + 8*yh^2) - 30*xh^2*(1 - 7*yh^2 + 8*yh^4) + yh^2*(10 - 33*yh^2 + 24*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[1001/(2*Pi)]*(-30*xh^2*yh^2*(1 - 6*zh^2 + 2*yh^2*(-1 + 8*zh^2)) + yh^4*(5 - 30*zh^2 + 6*yh^2*(-1 + 8*zh^2)) + 5*xh^4*(1 - 6*zh^2 + 6*yh^2*(-1 + 8*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[1001/(2*Pi)]*yh*(15*xh^5*(-1 + 8*zh^2) - 10*xh^3*(-1 + 6*zh^2 + 3*yh^2*(-1 + 8*zh^2)) + xh*yh^2*(-10 + 60*zh^2 + 3*yh^2*(-1 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[1001/(2*Pi)]*(-30*xh^2*yh^2*(1 - 6*zh^2 + 2*yh^2*(-1 + 8*zh^2)) + yh^4*(5 - 30*zh^2 + 6*yh^2*(-1 + 8*zh^2)) + 5*xh^4*(1 - 6*zh^2 + 6*yh^2*(-1 + 8*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdz = (9*Sqrt[1001/(2*Pi)]*yh*(5*xh^4 - 10*xh^2*yh^2 + yh^4)*zh*(-3 + 8*zh^2))/(8*r^2)
        y=2.366619162231752E0_r8*(yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh)
        gx=2.366619162231752E0_r8*(-2*ir*xh*yh*(15*xh**4 + 10*yh**2 + 3*yh**4 - 10*xh**2*(1 + 3*yh**2))*zh)
        gy=2.366619162231752E0_r8*(-(ir*(-5*xh**4 + 30*xh**2*(1 + xh**2)*yh**2 - 5*(1 + 12*xh**2)*yh**4 + 6*yh**6)*zh))
        gz=2.366619162231752E0_r8*(-(ir*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*(-1 + 6*zh**2)))
        gxx=2.366619162231752E0_r8*(2*ir**2*yh*(120*xh**6 - yh**2*(10 + 3*yh**2) - 15*xh**4*(9 + 16*yh**2) + 6*xh**2*(5 + 25*yh**2 + 4*yh**4))*zh)
        gxy=2.366619162231752E0_r8*(2*ir**2*(15*xh**5*(-1 + 8*yh**2) + xh**3*(10 + 30*yh**2 - 240*yh**4) + 3*xh*yh**2*(-10 + 15*yh**2 + 8*yh**4))*zh)
        gxz=2.366619162231752E0_r8*(2*ir**2*yh*(15*xh**5*(-1 + 8*zh**2) - 10*xh**3*(-1 + 6*zh**2 + 3*yh**2*(-1 + 8*zh**2)) + xh*yh**2*(-10 + 60*zh**2 + 3*yh**2*(-1 + 8*zh**2))))
        gyy=2.366619162231752E0_r8*(2*ir**2*yh*(15*xh**4*(-3 + 8*yh**2) - 30*xh**2*(1 - 7*yh**2 + 8*yh**4) + yh**2*(10 - 33*yh**2 + 24*yh**4))*zh)
        gyz=2.366619162231752E0_r8*(ir**2*(-30*xh**2*yh**2*(1 - 6*zh**2 + 2*yh**2*(-1 + 8*zh**2)) + yh**4*(5 - 30*zh**2 + 6*yh**2*(-1 + 8*zh**2)) + 5*xh**4*(1 - 6*zh**2 + 6*yh**2*(-1 + 8*zh**2))))
        gzz=2.366619162231752E0_r8*(6*ir**2*yh*(5*xh**4 - 10*xh**2*yh**2 + yh**4)*zh*(-3 + 8*zh**2))
    case(39)
        ! l=6, m=-4, iang=39 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*xh*yh*(xh^2 - yh^2)*(-1 + 11*zh^2))/8
        ! dYlm/dx = (3*Sqrt[91/Pi]*yh*(-5*xh^4 + 6*xh^6 + yh^4 - 6*xh^2*yh^4 - 10*(6*xh^4 + yh^2 - 3*xh^2*(1 + 2*yh^2))*zh^2))/(8*r)
        ! dYlm/dy = (3*Sqrt[91/Pi]*xh*(-xh^4 + 6*xh^4*yh^2 + 5*yh^4 - 6*yh^6 + 10*(xh^2 - 3*(1 + 2*xh^2)*yh^2 + 6*yh^4)*zh^2))/(8*r)
        ! dYlm/dz = (3*Sqrt[91/Pi]*xh*(xh - yh)*yh*(xh + yh)*zh*(10 + 3*xh^2 + 3*yh^2 - 30*zh^2))/(4*r)
        ! d^2Ylm/dxdx = (3*Sqrt[91/Pi]*xh*yh*(-24*xh^6 - 9*yh^4 + 30*zh^2 + 90*yh^2*zh^2 + 3*xh^4*(11 + 80*zh^2) + 2*xh^2*(-5 + 12*yh^4 - 105*zh^2 - 120*yh^2*zh^2)))/(4*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[91/Pi]*(xh^2 - yh^2)*(-6*yh^4 + 6*xh^4*(-1 + 8*yh^2) - 30*zh^2 + yh^2*(5 + 60*zh^2) + xh^2*(5 + 48*yh^4 + 60*zh^2 - 12*yh^2*(3 + 40*zh^2))))/(8*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[91/Pi]*yh*zh*(-24*xh^6 + 15*xh^4*(-3 + 16*zh^2) + yh^2*(-10 - 3*yh^2 + 30*zh^2) + 6*xh^2*(5 + 4*yh^4 - 15*zh^2 + yh^2*(10 - 40*zh^2))))/(4*r^2)
        ! d^2Ylm/dydx = (-3*Sqrt[91/Pi]*(xh^2 - yh^2)*(-6*yh^4 + 6*xh^4*(-1 + 8*yh^2) - 30*zh^2 + yh^2*(5 + 60*zh^2) + xh^2*(5 + 48*yh^4 + 60*zh^2 - 12*yh^2*(3 + 40*zh^2))))/(8*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[91/Pi]*xh*yh*(-24*yh^6 + 3*xh^4*(-3 + 8*yh^2) + 30*zh^2 - 30*xh^2*(-3 + 8*yh^2)*zh^2 - 10*yh^2*(1 + 21*zh^2) + 3*yh^4*(11 + 80*zh^2)))/(4*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[91/Pi]*xh*zh*(3*xh^4*(-1 + 8*yh^2) - 3*yh^2*(-10 + 8*yh^4 + 30*zh^2 + yh^2*(15 - 80*zh^2)) - 10*xh^2*(1 - 3*zh^2 + 6*yh^2*(-1 + 4*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[91/Pi]*yh*zh*(-24*xh^6 + 15*xh^4*(-3 + 16*zh^2) + yh^2*(-10 - 3*yh^2 + 30*zh^2) + 6*xh^2*(5 + 4*yh^4 - 15*zh^2 + yh^2*(10 - 40*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[91/Pi]*xh*zh*(3*xh^4*(-1 + 8*yh^2) - 3*yh^2*(-10 + 8*yh^4 + 30*zh^2 + yh^2*(15 - 80*zh^2)) - 10*xh^2*(1 - 3*zh^2 + 6*yh^2*(-1 + 4*zh^2))))/(4*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[91/Pi]*xh*yh*(xh^2 - yh^2)*(3*xh^2*(-1 + 8*zh^2) + 3*yh^2*(-1 + 8*zh^2) - 10*(1 - 15*zh^2 + 24*zh^4)))/(4*r^2)
        y=2.018259602914897E0_r8*(xh*yh*(xh**2 - yh**2)*(-1 + 11*zh**2))
        gx=2.018259602914897E0_r8*(ir*yh*(-5*xh**4 + 6*xh**6 + yh**4 - 6*xh**2*yh**4 - 10*(6*xh**4 + yh**2 - 3*xh**2*(1 + 2*yh**2))*zh**2))
        gy=2.018259602914897E0_r8*(ir*xh*(-xh**4 + 6*xh**4*yh**2 + 5*yh**4 - 6*yh**6 + 10*(xh**2 - 3*(1 + 2*xh**2)*yh**2 + 6*yh**4)*zh**2))
        gz=2.018259602914897E0_r8*(2*ir*xh*(xh - yh)*yh*(xh + yh)*zh*(10 + 3*xh**2 + 3*yh**2 - 30*zh**2))
        gxx=2.018259602914897E0_r8*(2*ir**2*xh*yh*(-24*xh**6 - 9*yh**4 + 30*zh**2 + 90*yh**2*zh**2 + 3*xh**4*(11 + 80*zh**2) + 2*xh**2*(-5 + 12*yh**4 - 105*zh**2 - 120*yh**2*zh**2)))
        gxy=2.018259602914897E0_r8*(-(ir**2*(xh**2 - yh**2)*(-6*yh**4 + 6*xh**4*(-1 + 8*yh**2) - 30*zh**2 + yh**2*(5 + 60*zh**2) + xh**2*(5 + 48*yh**4 + 60*zh**2 - 12*yh**2*(3 + 40*zh**2)))))
        gxz=2.018259602914897E0_r8*(2*ir**2*yh*zh*(-24*xh**6 + 15*xh**4*(-3 + 16*zh**2) + yh**2*(-10 - 3*yh**2 + 30*zh**2) + 6*xh**2*(5 + 4*yh**4 - 15*zh**2 + yh**2*(10 - 40*zh**2))))
        gyy=2.018259602914897E0_r8*(-2*ir**2*xh*yh*(-24*yh**6 + 3*xh**4*(-3 + 8*yh**2) + 30*zh**2 - 30*xh**2*(-3 + 8*yh**2)*zh**2 - 10*yh**2*(1 + 21*zh**2) + 3*yh**4*(11 + 80*zh**2)))
        gyz=2.018259602914897E0_r8*(-2*ir**2*xh*zh*(3*xh**4*(-1 + 8*yh**2) - 3*yh**2*(-10 + 8*yh**4 + 30*zh**2 + yh**2*(15 - 80*zh**2)) - 10*xh**2*(1 - 3*zh**2 + 6*yh**2*(-1 + 4*zh**2))))
        gzz=2.018259602914897E0_r8*(-2*ir**2*xh*yh*(xh**2 - yh**2)*(3*xh**2*(-1 + 8*zh**2) + 3*yh**2*(-1 + 8*zh**2) - 10*(1 - 15*zh**2 + 24*zh**4)))
    case(40)
        ! l=6, m=-3, iang=40 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(-3 + 11*zh^2))/16
        ! dYlm/dx = (-3*Sqrt[1365/(2*Pi)]*xh*yh*zh*(3 - 6*xh^2 + 2*yh^2 - 11*(1 - 3*xh^2 + yh^2)*zh^2))/(8*r)
        ! dYlm/dy = (3*Sqrt[1365/(2*Pi)]*zh*(xh^2*(-3 + 11*zh^2 + yh^2*(12 - 66*zh^2)) + yh^2*(3 - 11*zh^2 + yh^2*(-4 + 22*zh^2))))/(16*r)
        ! dYlm/dz = (3*Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*(1 - 15*zh^2 + 22*zh^4))/(16*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[1365/(2*Pi)]*yh*zh*(3 - 11*zh^2 + xh^4*(36 - 264*zh^2) + yh^2*(2 - 11*zh^2) + xh^2*(-30 + 165*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[1365/(2*Pi)]*xh*zh*(-3 + 11*zh^2 + yh^4*(12 - 88*zh^2) + yh^2*(6 - 33*zh^2) + 3*xh^2*(2 - 11*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[1365/(2*Pi)]*xh*yh*(-3 + 45*zh^2 - 66*zh^4 + yh^2*(-2 + 45*zh^2 - 88*zh^4) + 3*xh^2*(2 - 45*zh^2 + 88*zh^4)))/(8*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[1365/(2*Pi)]*xh*zh*(-3 + 11*zh^2 + yh^4*(12 - 88*zh^2) + yh^2*(6 - 33*zh^2) + 3*xh^2*(2 - 11*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[1365/(2*Pi)]*yh*zh*(-3 + 11*zh^2 + yh^2*(14 - 77*zh^2) + 4*yh^4*(-3 + 22*zh^2) - 3*xh^2*(6 - 33*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[1365/(2*Pi)]*(yh^2*(3 - 45*zh^2 + 66*zh^4 - 2*yh^2*(2 - 45*zh^2 + 88*zh^4)) + 3*xh^2*(-1 + 15*zh^2 - 22*zh^4 + 2*yh^2*(2 - 45*zh^2 + 88*zh^4))))/(16*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[1365/(2*Pi)]*xh*yh*(-3 + 45*zh^2 - 66*zh^4 + yh^2*(-2 + 45*zh^2 - 88*zh^4) + 3*xh^2*(2 - 45*zh^2 + 88*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[1365/(2*Pi)]*(yh^2*(3 - 45*zh^2 + 66*zh^4 - 2*yh^2*(2 - 45*zh^2 + 88*zh^4)) + 3*xh^2*(-1 + 15*zh^2 - 22*zh^4 + 2*yh^2*(2 - 45*zh^2 + 88*zh^4))))/(16*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[1365/(2*Pi)]*yh*(-3*xh^2 + yh^2)*zh*(17 - 89*zh^2 + 88*zh^4))/(8*r^2)
        y=-9.212052595149235E-1_r8*(yh*(-3*xh**2 + yh**2)*zh*(-3 + 11*zh**2))
        gx=-9.212052595149235E-1_r8*(6*ir*xh*yh*zh*(3 - 6*xh**2 + 2*yh**2 - 11*(1 - 3*xh**2 + yh**2)*zh**2))
        gy=-9.212052595149235E-1_r8*(-3*ir*zh*(xh**2*(-3 + 11*zh**2 + yh**2*(12 - 66*zh**2)) + yh**2*(3 - 11*zh**2 + yh**2*(-4 + 22*zh**2))))
        gz=-9.212052595149235E-1_r8*(-3*ir*yh*(-3*xh**2 + yh**2)*(1 - 15*zh**2 + 22*zh**4))
        gxx=-9.212052595149235E-1_r8*(6*ir**2*yh*zh*(3 - 11*zh**2 + xh**4*(36 - 264*zh**2) + yh**2*(2 - 11*zh**2) + xh**2*(-30 + 165*zh**2 + 4*yh**2*(-3 + 22*zh**2))))
        gxy=-9.212052595149235E-1_r8*(-6*ir**2*xh*zh*(-3 + 11*zh**2 + yh**4*(12 - 88*zh**2) + yh**2*(6 - 33*zh**2) + 3*xh**2*(2 - 11*zh**2 + 4*yh**2*(-3 + 22*zh**2))))
        gxz=-9.212052595149235E-1_r8*(-6*ir**2*xh*yh*(-3 + 45*zh**2 - 66*zh**4 + yh**2*(-2 + 45*zh**2 - 88*zh**4) + 3*xh**2*(2 - 45*zh**2 + 88*zh**4)))
        gyy=-9.212052595149235E-1_r8*(6*ir**2*yh*zh*(-3 + 11*zh**2 + yh**2*(14 - 77*zh**2) + 4*yh**4*(-3 + 22*zh**2) - 3*xh**2*(6 - 33*zh**2 + 4*yh**2*(-3 + 22*zh**2))))
        gyz=-9.212052595149235E-1_r8*(-3*ir**2*(yh**2*(3 - 45*zh**2 + 66*zh**4 - 2*yh**2*(2 - 45*zh**2 + 88*zh**4)) + 3*xh**2*(-1 + 15*zh**2 - 22*zh**4 + 2*yh**2*(2 - 45*zh**2 + 88*zh**4))))
        gzz=-9.212052595149235E-1_r8*(6*ir**2*yh*(-3*xh**2 + yh**2)*zh*(17 - 89*zh**2 + 88*zh**4))
    case(41)
        ! l=6, m=-2, iang=41 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*xh*yh*(1 - 18*zh^2 + 33*zh^4))/16
        ! dYlm/dx = (Sqrt[1365/(2*Pi)]*yh*(1 - 18*zh^2 + 33*zh^4 - 2*xh^2*(1 - 36*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dy = (Sqrt[1365/(2*Pi)]*xh*(1 - 18*zh^2 + 33*zh^4 - 2*yh^2*(1 - 36*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dz = -(Sqrt[1365/(2*Pi)]*xh*yh*zh*(19 - 102*zh^2 + 99*zh^4))/(8*r)
        ! d^2Ylm/dxdx = (Sqrt[1365/(2*Pi)]*yh*(-3*xh*(1 - 36*zh^2 + 99*zh^4) + 4*xh^3*(1 - 54*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dxdy = (Sqrt[1365/(2*Pi)]*(1 - 18*zh^2 + 33*zh^4 - 2*yh^2*(1 - 36*zh^2 + 99*zh^4) + 2*xh^2*(-1 + 36*zh^2 - 99*zh^4 + 4*yh^2*(1 - 54*zh^2 + 198*zh^4))))/(16*r^2)
        ! d^2Ylm/dxdz = (Sqrt[1365/(2*Pi)]*yh*zh*(-19 + 102*zh^2 - 99*zh^4 + xh^2*(76 - 612*zh^2 + 792*zh^4)))/(8*r^2)
        ! d^2Ylm/dydx = (Sqrt[1365/(2*Pi)]*(1 - 18*zh^2 + 33*zh^4 - 2*yh^2*(1 - 36*zh^2 + 99*zh^4) + 2*xh^2*(-1 + 36*zh^2 - 99*zh^4 + 4*yh^2*(1 - 54*zh^2 + 198*zh^4))))/(16*r^2)
        ! d^2Ylm/dydy = (Sqrt[1365/(2*Pi)]*xh*(-3*yh*(1 - 36*zh^2 + 99*zh^4) + 4*yh^3*(1 - 54*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dydz = (Sqrt[1365/(2*Pi)]*xh*zh*(-19 + 102*zh^2 - 99*zh^4 + yh^2*(76 - 612*zh^2 + 792*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdx = (Sqrt[1365/(2*Pi)]*yh*zh*(-19 + 102*zh^2 - 99*zh^4 + xh^2*(76 - 612*zh^2 + 792*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdy = (Sqrt[1365/(2*Pi)]*xh*zh*(-19 + 102*zh^2 - 99*zh^4 + yh^2*(76 - 612*zh^2 + 792*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdz = (Sqrt[1365/(2*Pi)]*xh*yh*(-19 + 382*zh^2 - 1107*zh^4 + 792*zh^6))/(8*r^2)
        y=9.212052595149235E-1_r8*(xh*yh*(1 - 18*zh**2 + 33*zh**4))
        gx=9.212052595149235E-1_r8*(ir*yh*(1 - 18*zh**2 + 33*zh**4 - 2*xh**2*(1 - 36*zh**2 + 99*zh**4)))
        gy=9.212052595149235E-1_r8*(ir*xh*(1 - 18*zh**2 + 33*zh**4 - 2*yh**2*(1 - 36*zh**2 + 99*zh**4)))
        gz=9.212052595149235E-1_r8*(-2*ir*xh*yh*zh*(19 - 102*zh**2 + 99*zh**4))
        gxx=9.212052595149235E-1_r8*(2*ir**2*yh*(-3*xh*(1 - 36*zh**2 + 99*zh**4) + 4*xh**3*(1 - 54*zh**2 + 198*zh**4)))
        gxy=9.212052595149235E-1_r8*(ir**2*(1 - 18*zh**2 + 33*zh**4 - 2*yh**2*(1 - 36*zh**2 + 99*zh**4) + 2*xh**2*(-1 + 36*zh**2 - 99*zh**4 + 4*yh**2*(1 - 54*zh**2 + 198*zh**4))))
        gxz=9.212052595149235E-1_r8*(2*ir**2*yh*zh*(-19 + 102*zh**2 - 99*zh**4 + xh**2*(76 - 612*zh**2 + 792*zh**4)))
        gyy=9.212052595149235E-1_r8*(2*ir**2*xh*(-3*yh*(1 - 36*zh**2 + 99*zh**4) + 4*yh**3*(1 - 54*zh**2 + 198*zh**4)))
        gyz=9.212052595149235E-1_r8*(2*ir**2*xh*zh*(-19 + 102*zh**2 - 99*zh**4 + yh**2*(76 - 612*zh**2 + 792*zh**4)))
        gzz=9.212052595149235E-1_r8*(2*ir**2*xh*yh*(-19 + 382*zh**2 - 1107*zh**4 + 792*zh**6))
    case(42)
        ! l=6, m=-1, iang=42 (=l^2+l+1+m)
        ! Ylm= (Sqrt[273/Pi]*yh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        ! dYlm/dx = -(Sqrt[273/Pi]*xh*yh*zh*(5 - 60*zh^2 + 99*zh^4))/(8*r)
        ! dYlm/dy = (Sqrt[273/Pi]*zh*(5 - 30*zh^2 + 33*zh^4 - 2*yh^2*(5 - 60*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dz = (Sqrt[273/Pi]*yh*(5 - 100*zh^2 + 285*zh^4 - 198*zh^6))/(16*r)
        ! d^2Ylm/dxdx = (Sqrt[273/Pi]*yh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*xh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dxdy = (Sqrt[273/Pi]*xh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*yh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dxdz = (Sqrt[273/Pi]*xh*yh*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6))/(8*r^2)
        ! d^2Ylm/dydx = (Sqrt[273/Pi]*xh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*yh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dydy = (Sqrt[273/Pi]*zh*(-3*yh*(5 - 60*zh^2 + 99*zh^4) + 4*yh^3*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dydz = (Sqrt[273/Pi]*(5 - 100*zh^2 + 285*zh^4 - 198*zh^6 + 2*yh^2*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6)))/(16*r^2)
        ! d^2Ylm/dzdx = (Sqrt[273/Pi]*xh*yh*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6))/(8*r^2)
        ! d^2Ylm/dzdy = (Sqrt[273/Pi]*(5 - 100*zh^2 + 285*zh^4 - 198*zh^6 + 2*yh^2*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6)))/(16*r^2)
        ! d^2Ylm/dzdz = (Sqrt[273/Pi]*yh*zh*(-105 + 770*zh^2 - 1449*zh^4 + 792*zh^6))/(8*r^2)
        y=5.826213625187314E-1_r8*(yh*zh*(5 - 30*zh**2 + 33*zh**4))
        gx=5.826213625187314E-1_r8*(-2*ir*xh*yh*zh*(5 - 60*zh**2 + 99*zh**4))
        gy=5.826213625187314E-1_r8*(ir*zh*(5 - 30*zh**2 + 33*zh**4 - 2*yh**2*(5 - 60*zh**2 + 99*zh**4)))
        gz=5.826213625187314E-1_r8*(ir*yh*(5 - 100*zh**2 + 285*zh**4 - 198*zh**6))
        gxx=5.826213625187314E-1_r8*(2*ir**2*yh*zh*(-5 + 60*zh**2 - 99*zh**4 + 4*xh**2*(5 - 90*zh**2 + 198*zh**4)))
        gxy=5.826213625187314E-1_r8*(2*ir**2*xh*zh*(-5 + 60*zh**2 - 99*zh**4 + 4*yh**2*(5 - 90*zh**2 + 198*zh**4)))
        gxz=5.826213625187314E-1_r8*(2*ir**2*xh*yh*(-5 + 200*zh**2 - 855*zh**4 + 792*zh**6))
        gyy=5.826213625187314E-1_r8*(2*ir**2*zh*(-3*yh*(5 - 60*zh**2 + 99*zh**4) + 4*yh**3*(5 - 90*zh**2 + 198*zh**4)))
        gyz=5.826213625187314E-1_r8*(ir**2*(5 - 100*zh**2 + 285*zh**4 - 198*zh**6 + 2*yh**2*(-5 + 200*zh**2 - 855*zh**4 + 792*zh**6)))
        gzz=5.826213625187314E-1_r8*(2*ir**2*yh*zh*(-105 + 770*zh**2 - 1449*zh**4 + 792*zh**6))
    case(43)
        ! l=6, m=0, iang=43 (=l^2+l+1+m)
        ! Ylm= (Sqrt[13/Pi]*(-5 + 105*zh^2 - 315*zh^4 + 231*zh^6))/32
        ! dYlm/dx = (-21*Sqrt[13/Pi]*xh*zh^2*(5 - 30*zh^2 + 33*zh^4))/(16*r)
        ! dYlm/dy = (-21*Sqrt[13/Pi]*yh*zh^2*(5 - 30*zh^2 + 33*zh^4))/(16*r)
        ! dYlm/dz = (-21*Sqrt[13/Pi]*zh*(-5 + 35*zh^2 - 63*zh^4 + 33*zh^6))/(16*r)
        ! d^2Ylm/dxdx = (21*Sqrt[13/Pi]*zh^2*(-5 + 30*zh^2 - 33*zh^4 + 4*xh^2*(5 - 45*zh^2 + 66*zh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (21*Sqrt[13/Pi]*xh*yh*zh^2*(5 - 45*zh^2 + 66*zh^4))/(4*r^2)
        ! d^2Ylm/dxdz = (21*Sqrt[13/Pi]*xh*zh*(-5 + 70*zh^2 - 189*zh^4 + 132*zh^6))/(8*r^2)
        ! d^2Ylm/dydx = (21*Sqrt[13/Pi]*xh*yh*zh^2*(5 - 45*zh^2 + 66*zh^4))/(4*r^2)
        ! d^2Ylm/dydy = (21*Sqrt[13/Pi]*zh^2*(-5 + 30*zh^2 - 33*zh^4 + 4*yh^2*(5 - 45*zh^2 + 66*zh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (21*Sqrt[13/Pi]*yh*zh*(-5 + 70*zh^2 - 189*zh^4 + 132*zh^6))/(8*r^2)
        ! d^2Ylm/dzdx = (21*Sqrt[13/Pi]*xh*zh*(-5 + 70*zh^2 - 189*zh^4 + 132*zh^6))/(8*r^2)
        ! d^2Ylm/dzdy = (21*Sqrt[13/Pi]*yh*zh*(-5 + 70*zh^2 - 189*zh^4 + 132*zh^6))/(8*r^2)
        ! d^2Ylm/dzdz = (21*Sqrt[13/Pi]*(5 - 115*zh^2 + 455*zh^4 - 609*zh^6 + 264*zh^8))/(16*r^2)
        y=6.356920226762843E-2_r8*(-5 + 105*zh**2 - 315*zh**4 + 231*zh**6)
        gx=6.356920226762843E-2_r8*(-42*ir*xh*zh**2*(5 - 30*zh**2 + 33*zh**4))
        gy=6.356920226762843E-2_r8*(-42*ir*yh*zh**2*(5 - 30*zh**2 + 33*zh**4))
        gz=6.356920226762843E-2_r8*(-42*ir*zh*(-5 + 35*zh**2 - 63*zh**4 + 33*zh**6))
        gxx=6.356920226762843E-2_r8*(42*ir**2*zh**2*(-5 + 30*zh**2 - 33*zh**4 + 4*xh**2*(5 - 45*zh**2 + 66*zh**4)))
        gxy=6.356920226762843E-2_r8*(168*ir**2*xh*yh*zh**2*(5 - 45*zh**2 + 66*zh**4))
        gxz=6.356920226762843E-2_r8*(84*ir**2*xh*zh*(-5 + 70*zh**2 - 189*zh**4 + 132*zh**6))
        gyy=6.356920226762843E-2_r8*(42*ir**2*zh**2*(-5 + 30*zh**2 - 33*zh**4 + 4*yh**2*(5 - 45*zh**2 + 66*zh**4)))
        gyz=6.356920226762843E-2_r8*(84*ir**2*yh*zh*(-5 + 70*zh**2 - 189*zh**4 + 132*zh**6))
        gzz=6.356920226762843E-2_r8*(42*ir**2*(5 - 115*zh**2 + 455*zh**4 - 609*zh**6 + 264*zh**8))
    case(44)
        ! l=6, m=1, iang=44 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[273/Pi]*xh*zh*(5 - 30*zh^2 + 33*zh^4))/16
        ! dYlm/dx = (Sqrt[273/Pi]*zh*(-5 + 30*zh^2 - 33*zh^4 + 2*xh^2*(5 - 60*zh^2 + 99*zh^4)))/(16*r)
        ! dYlm/dy = (Sqrt[273/Pi]*xh*yh*zh*(5 - 60*zh^2 + 99*zh^4))/(8*r)
        ! dYlm/dz = (Sqrt[273/Pi]*xh*(-5 + 100*zh^2 - 285*zh^4 + 198*zh^6))/(16*r)
        ! d^2Ylm/dxdx = -(Sqrt[273/Pi]*zh*(-3*xh*(5 - 60*zh^2 + 99*zh^4) + 4*xh^3*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dxdy = -(Sqrt[273/Pi]*yh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*xh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dxdz = (Sqrt[273/Pi]*(-5 + 100*zh^2 - 285*zh^4 + 198*zh^6 - 2*xh^2*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6)))/(16*r^2)
        ! d^2Ylm/dydx = -(Sqrt[273/Pi]*yh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*xh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dydy = -(Sqrt[273/Pi]*xh*zh*(-5 + 60*zh^2 - 99*zh^4 + 4*yh^2*(5 - 90*zh^2 + 198*zh^4)))/(8*r^2)
        ! d^2Ylm/dydz = -(Sqrt[273/Pi]*xh*yh*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6))/(8*r^2)
        ! d^2Ylm/dzdx = (Sqrt[273/Pi]*(-5 + 100*zh^2 - 285*zh^4 + 198*zh^6 - 2*xh^2*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6)))/(16*r^2)
        ! d^2Ylm/dzdy = -(Sqrt[273/Pi]*xh*yh*(-5 + 200*zh^2 - 855*zh^4 + 792*zh^6))/(8*r^2)
        ! d^2Ylm/dzdz = -(Sqrt[273/Pi]*xh*zh*(-105 + 770*zh^2 - 1449*zh^4 + 792*zh^6))/(8*r^2)
        y=-5.826213625187314E-1_r8*(xh*zh*(5 - 30*zh**2 + 33*zh**4))
        gx=-5.826213625187314E-1_r8*(-(ir*zh*(-5 + 30*zh**2 - 33*zh**4 + 2*xh**2*(5 - 60*zh**2 + 99*zh**4))))
        gy=-5.826213625187314E-1_r8*(-2*ir*xh*yh*zh*(5 - 60*zh**2 + 99*zh**4))
        gz=-5.826213625187314E-1_r8*(-(ir*xh*(-5 + 100*zh**2 - 285*zh**4 + 198*zh**6)))
        gxx=-5.826213625187314E-1_r8*(2*ir**2*zh*(-3*xh*(5 - 60*zh**2 + 99*zh**4) + 4*xh**3*(5 - 90*zh**2 + 198*zh**4)))
        gxy=-5.826213625187314E-1_r8*(2*ir**2*yh*zh*(-5 + 60*zh**2 - 99*zh**4 + 4*xh**2*(5 - 90*zh**2 + 198*zh**4)))
        gxz=-5.826213625187314E-1_r8*(-(ir**2*(-5 + 100*zh**2 - 285*zh**4 + 198*zh**6 - 2*xh**2*(-5 + 200*zh**2 - 855*zh**4 + 792*zh**6))))
        gyy=-5.826213625187314E-1_r8*(2*ir**2*xh*zh*(-5 + 60*zh**2 - 99*zh**4 + 4*yh**2*(5 - 90*zh**2 + 198*zh**4)))
        gyz=-5.826213625187314E-1_r8*(2*ir**2*xh*yh*(-5 + 200*zh**2 - 855*zh**4 + 792*zh**6))
        gzz=-5.826213625187314E-1_r8*(2*ir**2*xh*zh*(-105 + 770*zh**2 - 1449*zh**4 + 792*zh**6))
    case(45)
        ! l=6, m=2, iang=45 (=l^2+l+1+m)
        ! Ylm= (Sqrt[1365/(2*Pi)]*(xh^2 - yh^2)*(1 - 18*zh^2 + 33*zh^4))/32
        ! dYlm/dx = (Sqrt[1365/(2*Pi)]*xh*(1 - xh^2 + yh^2 + 18*(-1 + 2*xh^2 - 2*yh^2)*zh^2 + 33*(1 - 3*xh^2 + 3*yh^2)*zh^4))/(16*r)
        ! dYlm/dy = (Sqrt[1365/(2*Pi)]*yh*(-1 - xh^2 + yh^2 + 18*(1 + 2*xh^2 - 2*yh^2)*zh^2 - 33*(1 + 3*xh^2 - 3*yh^2)*zh^4))/(16*r)
        ! dYlm/dz = -(Sqrt[1365/(2*Pi)]*(xh - yh)*(xh + yh)*zh*(19 - 102*zh^2 + 99*zh^4))/(16*r)
        ! d^2Ylm/dxdx = (Sqrt[1365/(2*Pi)]*(1 - 18*zh^2 + 33*zh^4 + yh^2*(1 - 36*zh^2 + 99*zh^4) + 4*xh^4*(1 - 54*zh^2 + 198*zh^4) - xh^2*(5*(1 - 36*zh^2 + 99*zh^4) + 4*yh^2*(1 - 54*zh^2 + 198*zh^4))))/(16*r^2)
        ! d^2Ylm/dxdy = (Sqrt[1365/(2*Pi)]*xh*yh*(xh^2 - yh^2)*(1 - 54*zh^2 + 198*zh^4))/(4*r^2)
        ! d^2Ylm/dxdz = (Sqrt[1365/(2*Pi)]*xh*zh*(-19 + 102*zh^2 - 99*zh^4 + yh^2*(-38 + 306*zh^2 - 396*zh^4) + xh^2*(38 - 306*zh^2 + 396*zh^4)))/(8*r^2)
        ! d^2Ylm/dydx = (Sqrt[1365/(2*Pi)]*xh*yh*(xh^2 - yh^2)*(1 - 54*zh^2 + 198*zh^4))/(4*r^2)
        ! d^2Ylm/dydy = (Sqrt[1365/(2*Pi)]*(-1 + 18*zh^2 - 33*zh^4 + 5*yh^2*(1 - 36*zh^2 + 99*zh^4) - 4*yh^4*(1 - 54*zh^2 + 198*zh^4) + xh^2*(-1 + 36*zh^2 - 99*zh^4 + 4*yh^2*(1 - 54*zh^2 + 198*zh^4))))/(16*r^2)
        ! d^2Ylm/dydz = -(Sqrt[1365/(2*Pi)]*yh*zh*(-19 + 102*zh^2 - 99*zh^4 + xh^2*(-38 + 306*zh^2 - 396*zh^4) + yh^2*(38 - 306*zh^2 + 396*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdx = (Sqrt[1365/(2*Pi)]*xh*zh*(-19 + 102*zh^2 - 99*zh^4 + yh^2*(-38 + 306*zh^2 - 396*zh^4) + xh^2*(38 - 306*zh^2 + 396*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdy = -(Sqrt[1365/(2*Pi)]*yh*zh*(-19 + 102*zh^2 - 99*zh^4 + xh^2*(-38 + 306*zh^2 - 396*zh^4) + yh^2*(38 - 306*zh^2 + 396*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdz = (Sqrt[1365/(2*Pi)]*(xh^2 - yh^2)*(-19 + 382*zh^2 - 1107*zh^4 + 792*zh^6))/(16*r^2)
        y=4.606026297574617E-1_r8*((xh**2 - yh**2)*(1 - 18*zh**2 + 33*zh**4))
        gx=4.606026297574617E-1_r8*(2*ir*xh*(1 - xh**2 + yh**2 + 18*(-1 + 2*xh**2 - 2*yh**2)*zh**2 + 33*(1 - 3*xh**2 + 3*yh**2)*zh**4))
        gy=4.606026297574617E-1_r8*(2*ir*yh*(-1 - xh**2 + yh**2 + 18*(1 + 2*xh**2 - 2*yh**2)*zh**2 - 33*(1 + 3*xh**2 - 3*yh**2)*zh**4))
        gz=4.606026297574617E-1_r8*(-2*ir*(xh - yh)*(xh + yh)*zh*(19 - 102*zh**2 + 99*zh**4))
        gxx=4.606026297574617E-1_r8*(2*ir**2*(1 - 18*zh**2 + 33*zh**4 + yh**2*(1 - 36*zh**2 + 99*zh**4) + 4*xh**4*(1 - 54*zh**2 + 198*zh**4) - xh**2*(5*(1 - 36*zh**2 + 99*zh**4) + 4*yh**2*(1 - 54*zh**2 + 198*zh**4))))
        gxy=4.606026297574617E-1_r8*(8*ir**2*xh*yh*(xh**2 - yh**2)*(1 - 54*zh**2 + 198*zh**4))
        gxz=4.606026297574617E-1_r8*(4*ir**2*xh*zh*(-19 + 102*zh**2 - 99*zh**4 + yh**2*(-38 + 306*zh**2 - 396*zh**4) + xh**2*(38 - 306*zh**2 + 396*zh**4)))
        gyy=4.606026297574617E-1_r8*(2*ir**2*(-1 + 18*zh**2 - 33*zh**4 + 5*yh**2*(1 - 36*zh**2 + 99*zh**4) - 4*yh**4*(1 - 54*zh**2 + 198*zh**4) + xh**2*(-1 + 36*zh**2 - 99*zh**4 + 4*yh**2*(1 - 54*zh**2 + 198*zh**4))))
        gyz=4.606026297574617E-1_r8*(-4*ir**2*yh*zh*(-19 + 102*zh**2 - 99*zh**4 + xh**2*(-38 + 306*zh**2 - 396*zh**4) + yh**2*(38 - 306*zh**2 + 396*zh**4)))
        gzz=4.606026297574617E-1_r8*(2*ir**2*(xh**2 - yh**2)*(-19 + 382*zh**2 - 1107*zh**4 + 792*zh**6))
    case(46)
        ! l=6, m=3, iang=46 (=l^2+l+1+m)
        ! Ylm= -(Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(-3 + 11*zh^2))/16
        ! dYlm/dx = (3*Sqrt[1365/(2*Pi)]*zh*(yh^2*(-3 + 11*zh^2) + xh^4*(-4 + 22*zh^2) + xh^2*(3 + 12*yh^2 - 11*(1 + 6*yh^2)*zh^2)))/(16*r)
        ! dYlm/dy = (3*Sqrt[1365/(2*Pi)]*xh*yh*zh*(-3 - 2*xh^2 + 6*yh^2 + 11*(1 + xh^2 - 3*yh^2)*zh^2))/(8*r)
        ! dYlm/dz = (3*Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*(1 - 15*zh^2 + 22*zh^4))/(16*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[1365/(2*Pi)]*xh*zh*(-3 + 11*zh^2 + 9*yh^2*(-2 + 11*zh^2) + 4*xh^4*(-3 + 22*zh^2) + xh^2*(14 - 77*zh^2 + yh^2*(36 - 264*zh^2))))/(8*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[1365/(2*Pi)]*yh*zh*(-3 + 11*zh^2 + xh^4*(12 - 88*zh^2) + yh^2*(6 - 33*zh^2) + 3*xh^2*(2 - 11*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[1365/(2*Pi)]*(3*yh^2*(1 - 15*zh^2 + 22*zh^4) + 2*xh^4*(2 - 45*zh^2 + 88*zh^4) - 3*xh^2*(1 - 15*zh^2 + 22*zh^4 + 2*yh^2*(2 - 45*zh^2 + 88*zh^4))))/(16*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[1365/(2*Pi)]*yh*zh*(-3 + 11*zh^2 + xh^4*(12 - 88*zh^2) + yh^2*(6 - 33*zh^2) + 3*xh^2*(2 - 11*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[1365/(2*Pi)]*xh*zh*(3 - 11*zh^2 + yh^4*(36 - 264*zh^2) + 15*yh^2*(-2 + 11*zh^2) + xh^2*(2 - 11*zh^2 + 4*yh^2*(-3 + 22*zh^2))))/(8*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[1365/(2*Pi)]*xh*yh*(3 - 45*zh^2 + 66*zh^4 + xh^2*(2 - 45*zh^2 + 88*zh^4) - 3*yh^2*(2 - 45*zh^2 + 88*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdx = (-3*Sqrt[1365/(2*Pi)]*(3*yh^2*(1 - 15*zh^2 + 22*zh^4) + 2*xh^4*(2 - 45*zh^2 + 88*zh^4) - 3*xh^2*(1 - 15*zh^2 + 22*zh^4 + 2*yh^2*(2 - 45*zh^2 + 88*zh^4))))/(16*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[1365/(2*Pi)]*xh*yh*(3 - 45*zh^2 + 66*zh^4 + xh^2*(2 - 45*zh^2 + 88*zh^4) - 3*yh^2*(2 - 45*zh^2 + 88*zh^4)))/(8*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[1365/(2*Pi)]*xh*(xh^2 - 3*yh^2)*zh*(17 - 89*zh^2 + 88*zh^4))/(8*r^2)
        y=-9.212052595149235E-1_r8*(xh*(xh**2 - 3*yh**2)*zh*(-3 + 11*zh**2))
        gx=-9.212052595149235E-1_r8*(-3*ir*zh*(yh**2*(-3 + 11*zh**2) + xh**4*(-4 + 22*zh**2) + xh**2*(3 + 12*yh**2 - 11*(1 + 6*yh**2)*zh**2)))
        gy=-9.212052595149235E-1_r8*(-6*ir*xh*yh*zh*(-3 - 2*xh**2 + 6*yh**2 + 11*(1 + xh**2 - 3*yh**2)*zh**2))
        gz=-9.212052595149235E-1_r8*(-3*ir*xh*(xh**2 - 3*yh**2)*(1 - 15*zh**2 + 22*zh**4))
        gxx=-9.212052595149235E-1_r8*(6*ir**2*xh*zh*(-3 + 11*zh**2 + 9*yh**2*(-2 + 11*zh**2) + 4*xh**4*(-3 + 22*zh**2) + xh**2*(14 - 77*zh**2 + yh**2*(36 - 264*zh**2))))
        gxy=-9.212052595149235E-1_r8*(-6*ir**2*yh*zh*(-3 + 11*zh**2 + xh**4*(12 - 88*zh**2) + yh**2*(6 - 33*zh**2) + 3*xh**2*(2 - 11*zh**2 + 4*yh**2*(-3 + 22*zh**2))))
        gxz=-9.212052595149235E-1_r8*(3*ir**2*(3*yh**2*(1 - 15*zh**2 + 22*zh**4) + 2*xh**4*(2 - 45*zh**2 + 88*zh**4) - 3*xh**2*(1 - 15*zh**2 + 22*zh**4 + 2*yh**2*(2 - 45*zh**2 + 88*zh**4))))
        gyy=-9.212052595149235E-1_r8*(6*ir**2*xh*zh*(3 - 11*zh**2 + yh**4*(36 - 264*zh**2) + 15*yh**2*(-2 + 11*zh**2) + xh**2*(2 - 11*zh**2 + 4*yh**2*(-3 + 22*zh**2))))
        gyz=-9.212052595149235E-1_r8*(6*ir**2*xh*yh*(3 - 45*zh**2 + 66*zh**4 + xh**2*(2 - 45*zh**2 + 88*zh**4) - 3*yh**2*(2 - 45*zh**2 + 88*zh**4)))
        gzz=-9.212052595149235E-1_r8*(6*ir**2*xh*(xh**2 - 3*yh**2)*zh*(17 - 89*zh**2 + 88*zh**4))
    case(47)
        ! l=6, m=4, iang=47 (=l^2+l+1+m)
        ! Ylm= (3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(-1 + 11*zh^2))/32
        ! dYlm/dx = (3*Sqrt[91/Pi]*xh*(3*xh^6 + yh^4*(5 + 3*yh^2) - 30*yh^2*(2 + yh^2)*zh^2 - 3*xh^4*(1 + 5*yh^2 + 10*zh^2) + 5*xh^2*(-3*yh^4 + 4*zh^2 + yh^2*(2 + 36*zh^2))))/(16*r)
        ! dYlm/dy = (3*Sqrt[91/Pi]*yh*(3*xh^6 + 3*yh^6 + 20*yh^2*zh^2 - 5*xh^4*(-1 + 3*yh^2 + 6*zh^2) - 3*yh^4*(1 + 10*zh^2) + 5*xh^2*(-3*yh^4 - 12*zh^2 + yh^2*(2 + 36*zh^2))))/(16*r)
        ! dYlm/dz = (3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*zh*(10 + 3*xh^2 + 3*yh^2 - 30*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[91/Pi]*(24*xh^8 - 3*yh^6 + 60*yh^2*zh^2 + 5*yh^4*(-1 + 6*zh^2) - 3*xh^6*(13 + 40*yh^2 + 80*zh^2) &
        !               + 3*xh^2*(8*yh^6 - 20*zh^2 + yh^4*(25 - 80*zh^2) - 10*yh^2*(1 + 30*zh^2)) - 15*xh^4*(-1 + 8*yh^4 - 18*zh^2 - 3*yh^2*(3 + 32*zh^2))))/(16*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[91/Pi]*xh*yh*(6*xh^6 + 6*yh^6 + 30*zh^2 + yh^4*(3 - 60*zh^2) + xh^4*(3 - 30*yh^2 - 60*zh^2) - 5*yh^2*(1 + 12*zh^2) - 5*xh^2*(1 + 6*yh^4 + 12*zh^2 - 6*yh^2*(1 + 12*zh^2))))/(4*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[91/Pi]*zh*(-12*xh^7 + 3*xh^5*(-7 + 20*yh^2 + 40*zh^2) + 10*xh^3*(2 + 6*yh^4 - 6*zh^2 + yh^2*(15 - 72*zh^2)) - 3*xh*yh^2*(20 + 4*yh^4 - 60*zh^2 - 5*yh^2*(-3 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dydx = (-3*Sqrt[91/Pi]*xh*yh*(6*xh^6 + 6*yh^6 + 30*zh^2 + yh^4*(3 - 60*zh^2) + xh^4*(3 - 30*yh^2 - 60*zh^2) - 5*yh^2*(1 + 12*zh^2) - 5*xh^2*(1 + 6*yh^4 + 12*zh^2 - 6*yh^2*(1 + 12*zh^2))))/(4*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[91/Pi]*(3*xh^6*(-1 + 8*yh^2) - 5*xh^4*(1 + 24*yh^4 - 6*zh^2 + 3*yh^2*(-5 + 16*zh^2)) &
        !- 15*xh^2*(8*yh^6 - 4*zh^2 - 3*yh^4*(3 + 32*zh^2) + yh^2*(2 + 60*zh^2)) + 3*(8*yh^8 - 20*yh^2*zh^2 - yh^6*(13 + 80*zh^2) + yh^4*(5 + 90*zh^2))))/(16*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[91/Pi]*yh*zh*(-12*xh^6 + 15*xh^4*(-3 + 4*yh^2 + 8*zh^2) + 30*xh^2*(-2 + 2*yh^4 + 6*zh^2 + yh^2*(5 - 24*zh^2)) + yh^2*(20 - 12*yh^4 - 60*zh^2 + 3*yh^2*(-7 + 40*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[91/Pi]*zh*(-12*xh^7 + 3*xh^5*(-7 + 20*yh^2 + 40*zh^2) + 10*xh^3*(2 + 6*yh^4 - 6*zh^2 + yh^2*(15 - 72*zh^2)) - 3*xh*yh^2*(20 + 4*yh^4 - 60*zh^2 - 5*yh^2*(-3 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[91/Pi]*yh*zh*(-12*xh^6 + 15*xh^4*(-3 + 4*yh^2 + 8*zh^2) + 30*xh^2*(-2 + 2*yh^4 + 6*zh^2 + yh^2*(5 - 24*zh^2)) + yh^2*(20 - 12*yh^4 - 60*zh^2 + 3*yh^2*(-7 + 40*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdz = (-3*Sqrt[91/Pi]*(xh^4 - 6*xh^2*yh^2 + yh^4)*(3*xh^2*(-1 + 8*zh^2) + 3*yh^2*(-1 + 8*zh^2) - 10*(1 - 15*zh^2 + 24*zh^4)))/(16*r^2)
        y=5.045649007287242E-1_r8*((xh**4 - 6*xh**2*yh**2 + yh**4)*(-1 + 11*zh**2))
        gx=5.045649007287242E-1_r8*(2*ir*xh*(3*xh**6 + yh**4*(5 + 3*yh**2) - 30*yh**2*(2 + yh**2)*zh**2 - 3*xh**4*(1 + 5*yh**2 + 10*zh**2) + 5*xh**2*(-3*yh**4 + 4*zh**2 + yh**2*(2 + 36*zh**2))))
        gy=5.045649007287242E-1_r8*(2*ir*yh*(3*xh**6 + 3*yh**6 + 20*yh**2*zh**2 - 5*xh**4*(-1 + 3*yh**2 + 6*zh**2) - 3*yh**4*(1 + 10*zh**2) + 5*xh**2*(-3*yh**4 - 12*zh**2 + yh**2*(2 + 36*zh**2))))
        gz=5.045649007287242E-1_r8*(2*ir*(xh**4 - 6*xh**2*yh**2 + yh**4)*zh*(10 + 3*xh**2 + 3*yh**2 - 30*zh**2))
        gxx=5.045649007287242E-1_r8*(-2*ir**2*(24*xh**8 - 3*yh**6 + 60*yh**2*zh**2 + 5*yh**4*(-1 + 6*zh**2) - 3*xh**6*(13 + 40*yh**2 + 80*zh**2) &
            + 3*xh**2*(8*yh**6 - 20*zh**2 + yh**4*(25 - 80*zh**2) - 10*yh**2*(1 + 30*zh**2)) - 15*xh**4*(-1 + 8*yh**4 - 18*zh**2 - 3*yh**2*(3 + 32*zh**2))))
        gxy=5.045649007287242E-1_r8*(-8*ir**2*xh*yh*(6*xh**6 + 6*yh**6 + 30*zh**2 + yh**4*(3 - 60*zh**2) + xh**4*(3 - 30*yh**2 - 60*zh**2) - 5*yh**2*(1 + 12*zh**2) - 5*xh**2*(1 + 6*yh**4 + 12*zh**2 - 6*yh**2*(1 + 12*zh**2))))
        gxz=5.045649007287242E-1_r8*(4*ir**2*zh*(-12*xh**7 + 3*xh**5*(-7 + 20*yh**2 + 40*zh**2) + 10*xh**3*(2 + 6*yh**4 - 6*zh**2 + yh**2*(15 - 72*zh**2)) - 3*xh*yh**2*(20 + 4*yh**4 - 60*zh**2 - 5*yh**2*(-3 + 8*zh**2))))
        gyy=5.045649007287242E-1_r8*(-2*ir**2*(3*xh**6*(-1 + 8*yh**2) - 5*xh**4*(1 + 24*yh**4 - 6*zh**2 + 3*yh**2*(-5 + 16*zh**2)) &
           - 15*xh**2*(8*yh**6 - 4*zh**2 - 3*yh**4*(3 + 32*zh**2) + yh**2*(2 + 60*zh**2)) + 3*(8*yh**8 - 20*yh**2*zh**2 - yh**6*(13 + 80*zh**2) + yh**4*(5 + 90*zh**2))))
        gyz=5.045649007287242E-1_r8*(4*ir**2*yh*zh*(-12*xh**6 + 15*xh**4*(-3 + 4*yh**2 + 8*zh**2) + 30*xh**2*(-2 + 2*yh**4 + 6*zh**2 + yh**2*(5 - 24*zh**2)) + yh**2*(20 - 12*yh**4 - 60*zh**2 + 3*yh**2*(-7 + 40*zh**2))))
        gzz=5.045649007287242E-1_r8*(-2*ir**2*(xh**4 - 6*xh**2*yh**2 + yh**4)*(3*xh**2*(-1 + 8*zh**2) + 3*yh**2*(-1 + 8*zh**2) - 10*(1 - 15*zh**2 + 24*zh**4)))
    case(48)
        ! l=6, m=5, iang=48 (=l^2+l+1+m)
        ! Ylm= (-3*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh)/16
        ! dYlm/dx = (3*Sqrt[1001/(2*Pi)]*(6*xh^6 - 5*yh^4 - 5*xh^4*(1 + 12*yh^2) + 30*xh^2*(yh^2 + yh^4))*zh)/(16*r)
        ! dYlm/dy = (3*Sqrt[1001/(2*Pi)]*xh*yh*(3*xh^4 + xh^2*(10 - 30*yh^2) + 5*yh^2*(-2 + 3*yh^2))*zh)/(8*r)
        ! dYlm/dz = (3*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*(-1 + 6*zh^2))/(16*r)
        ! d^2Ylm/dxdx = (-3*Sqrt[1001/(2*Pi)]*(24*xh^7 - 15*xh*yh^2*(2 + 3*yh^2) - 3*xh^5*(11 + 80*yh^2) + 10*xh^3*(1 + 21*yh^2 + 12*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dxdy = (-3*Sqrt[1001/(2*Pi)]*yh*(24*xh^6 + xh^4*(45 - 240*yh^2) + 5*yh^2*(2 - 3*yh^2) + 30*xh^2*(-1 + yh^2 + 4*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dxdz = (-3*Sqrt[1001/(2*Pi)]*(5*yh^4*(1 - 6*zh^2) + 6*xh^6*(-1 + 8*zh^2) + 30*xh^2*yh^2*(-1 + 6*zh^2 + yh^2*(-1 + 8*zh^2)) - 5*xh^4*(-1 + 6*zh^2 + 12*yh^2*(-1 + 8*zh^2))))/(16*r^2)
        ! d^2Ylm/dydx = (-3*Sqrt[1001/(2*Pi)]*yh*(24*xh^6 + xh^4*(45 - 240*yh^2) + 5*yh^2*(2 - 3*yh^2) + 30*xh^2*(-1 + yh^2 + 4*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dydy = (-3*Sqrt[1001/(2*Pi)]*xh*(3*xh^4*(-1 + 8*yh^2) + 15*yh^2*(2 - 9*yh^2 + 8*yh^4) - 10*xh^2*(1 - 15*yh^2 + 24*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dydz = (-3*Sqrt[1001/(2*Pi)]*xh*yh*(3*xh^4*(-1 + 8*zh^2) + 5*yh^2*(2 - 12*zh^2 + 3*yh^2*(-1 + 8*zh^2)) - 10*xh^2*(1 - 6*zh^2 + 3*yh^2*(-1 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdx = (-3*Sqrt[1001/(2*Pi)]*(5*yh^4*(1 - 6*zh^2) + 6*xh^6*(-1 + 8*zh^2) + 30*xh^2*yh^2*(-1 + 6*zh^2 + yh^2*(-1 + 8*zh^2)) - 5*xh^4*(-1 + 6*zh^2 + 12*yh^2*(-1 + 8*zh^2))))/(16*r^2)
        ! d^2Ylm/dzdy = (-3*Sqrt[1001/(2*Pi)]*xh*yh*(3*xh^4*(-1 + 8*zh^2) + 5*yh^2*(2 - 12*zh^2 + 3*yh^2*(-1 + 8*zh^2)) - 10*xh^2*(1 - 6*zh^2 + 3*yh^2*(-1 + 8*zh^2))))/(8*r^2)
        ! d^2Ylm/dzdz = (-9*Sqrt[1001/(2*Pi)]*xh*(xh^4 - 10*xh^2*yh^2 + 5*yh^4)*zh*(-3 + 8*zh^2))/(8*r^2)
        y=-2.366619162231752E0_r8*(xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh)
        gx=-2.366619162231752E0_r8*(-(ir*(6*xh**6 - 5*yh**4 - 5*xh**4*(1 + 12*yh**2) + 30*xh**2*(yh**2 + yh**4))*zh))
        gy=-2.366619162231752E0_r8*(-2*ir*xh*yh*(3*xh**4 + xh**2*(10 - 30*yh**2) + 5*yh**2*(-2 + 3*yh**2))*zh)
        gz=-2.366619162231752E0_r8*(-(ir*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*(-1 + 6*zh**2)))
        gxx=-2.366619162231752E0_r8*(2*ir**2*(24*xh**7 - 15*xh*yh**2*(2 + 3*yh**2) - 3*xh**5*(11 + 80*yh**2) + 10*xh**3*(1 + 21*yh**2 + 12*yh**4))*zh)
        gxy=-2.366619162231752E0_r8*(2*ir**2*yh*(24*xh**6 + xh**4*(45 - 240*yh**2) + 5*yh**2*(2 - 3*yh**2) + 30*xh**2*(-1 + yh**2 + 4*yh**4))*zh)
        gxz=-2.366619162231752E0_r8*(ir**2*(5*yh**4*(1 - 6*zh**2) + 6*xh**6*(-1 + 8*zh**2) + 30*xh**2*yh**2*(-1 + 6*zh**2 + yh**2*(-1 + 8*zh**2)) - 5*xh**4*(-1 + 6*zh**2 + 12*yh**2*(-1 + 8*zh**2))))
        gyy=-2.366619162231752E0_r8*(2*ir**2*xh*(3*xh**4*(-1 + 8*yh**2) + 15*yh**2*(2 - 9*yh**2 + 8*yh**4) - 10*xh**2*(1 - 15*yh**2 + 24*yh**4))*zh)
        gyz=-2.366619162231752E0_r8*(2*ir**2*xh*yh*(3*xh**4*(-1 + 8*zh**2) + 5*yh**2*(2 - 12*zh**2 + 3*yh**2*(-1 + 8*zh**2)) - 10*xh**2*(1 - 6*zh**2 + 3*yh**2*(-1 + 8*zh**2))))
        gzz=-2.366619162231752E0_r8*(6*ir**2*xh*(xh**4 - 10*xh**2*yh**2 + 5*yh**4)*zh*(-3 + 8*zh**2))
    case(49)
        ! l=6, m=6, iang=49 (=l^2+l+1+m)
        ! Ylm= (Sqrt[3003/(2*Pi)]*(xh^6 - 15*xh^4*yh^2 + 15*xh^2*yh^4 - yh^6))/32
        ! dYlm/dx = (-3*Sqrt[3003/(2*Pi)]*(xh^7 - xh*yh^4*(5 + yh^2) + 5*xh^3*yh^2*(2 + 3*yh^2) - xh^5*(1 + 15*yh^2)))/(16*r)
        ! dYlm/dy = (-3*Sqrt[3003/(2*Pi)]*yh*(xh^6 + yh^4 - yh^6 + xh^4*(5 - 15*yh^2) + 5*xh^2*yh^2*(-2 + 3*yh^2)))/(16*r)
        ! dYlm/dz = (3*Sqrt[3003/(2*Pi)]*(-xh^6 + 15*xh^4*yh^2 - 15*xh^2*yh^4 + yh^6)*zh)/(16*r)
        ! d^2Ylm/dxdx = (3*Sqrt[3003/(2*Pi)]*(8*xh^8 + yh^4*(5 + yh^2) - xh^6*(13 + 120*yh^2) - xh^2*yh^2*(30 + 75*yh^2 + 8*yh^4) + 5*xh^4*(1 + 27*yh^2 + 24*yh^4)))/(16*r^2)
        ! d^2Ylm/dxdy = (3*Sqrt[3003/(2*Pi)]*xh*yh*(2*xh^6 + 5*yh^2 - 6*yh^4 - 2*yh^6 + xh^4*(6 - 30*yh^2) + 5*xh^2*(-1 + 6*yh^4)))/(4*r^2)
        ! d^2Ylm/dxdz = (3*Sqrt[3003/(2*Pi)]*(4*xh^7 - xh*yh^4*(15 + 4*yh^2) - 3*xh^5*(1 + 20*yh^2) + 30*xh^3*(yh^2 + 2*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dydx = (3*Sqrt[3003/(2*Pi)]*xh*yh*(2*xh^6 + 5*yh^2 - 6*yh^4 - 2*yh^6 + xh^4*(6 - 30*yh^2) + 5*xh^2*(-1 + 6*yh^4)))/(4*r^2)
        ! d^2Ylm/dydy = (3*Sqrt[3003/(2*Pi)]*(-5*yh^4 + 13*yh^6 - 8*yh^8 + xh^6*(-1 + 8*yh^2) + 15*xh^2*yh^2*(2 - 9*yh^2 + 8*yh^4) - 5*xh^4*(1 - 15*yh^2 + 24*yh^4)))/(16*r^2)
        ! d^2Ylm/dydz = (3*Sqrt[3003/(2*Pi)]*yh*(4*xh^6 + 3*yh^4 - 4*yh^6 + xh^4*(15 - 60*yh^2) + 30*xh^2*yh^2*(-1 + 2*yh^2))*zh)/(8*r^2)
        ! d^2Ylm/dzdx = (3*Sqrt[3003/(2*Pi)]*(4*xh^7 - xh*yh^4*(15 + 4*yh^2) - 3*xh^5*(1 + 20*yh^2) + 30*xh^3*(yh^2 + 2*yh^4))*zh)/(8*r^2)
        ! d^2Ylm/dzdy = (3*Sqrt[3003/(2*Pi)]*yh*(4*xh^6 + 3*yh^4 - 4*yh^6 + xh^4*(15 - 60*yh^2) + 30*xh^2*yh^2*(-1 + 2*yh^2))*zh)/(8*r^2)
        ! d^2Ylm/dzdz = (3*Sqrt[3003/(2*Pi)]*(xh^6 - 15*xh^4*yh^2 + 15*xh^2*yh^4 - yh^6)*(-1 + 8*zh^2))/(16*r^2)
        y=6.831841051919143E-1_r8*(xh**6 - 15*xh**4*yh**2 + 15*xh**2*yh**4 - yh**6)
        gx=6.831841051919143E-1_r8*(-6*ir*(xh**7 - xh*yh**4*(5 + yh**2) + 5*xh**3*yh**2*(2 + 3*yh**2) - xh**5*(1 + 15*yh**2)))
        gy=6.831841051919143E-1_r8*(-6*ir*yh*(xh**6 + yh**4 - yh**6 + xh**4*(5 - 15*yh**2) + 5*xh**2*yh**2*(-2 + 3*yh**2)))
        gz=6.831841051919143E-1_r8*(6*ir*(-xh**6 + 15*xh**4*yh**2 - 15*xh**2*yh**4 + yh**6)*zh)
        gxx=6.831841051919143E-1_r8*(6*ir**2*(8*xh**8 + yh**4*(5 + yh**2) - xh**6*(13 + 120*yh**2) - xh**2*yh**2*(30 + 75*yh**2 + 8*yh**4) + 5*xh**4*(1 + 27*yh**2 + 24*yh**4)))
        gxy=6.831841051919143E-1_r8*(24*ir**2*xh*yh*(2*xh**6 + 5*yh**2 - 6*yh**4 - 2*yh**6 + xh**4*(6 - 30*yh**2) + 5*xh**2*(-1 + 6*yh**4)))
        gxz=6.831841051919143E-1_r8*(12*ir**2*(4*xh**7 - xh*yh**4*(15 + 4*yh**2) - 3*xh**5*(1 + 20*yh**2) + 30*xh**3*(yh**2 + 2*yh**4))*zh)
        gyy=6.831841051919143E-1_r8*(6*ir**2*(-5*yh**4 + 13*yh**6 - 8*yh**8 + xh**6*(-1 + 8*yh**2) + 15*xh**2*yh**2*(2 - 9*yh**2 + 8*yh**4) - 5*xh**4*(1 - 15*yh**2 + 24*yh**4)))
        gyz=6.831841051919143E-1_r8*(12*ir**2*yh*(4*xh**6 + 3*yh**4 - 4*yh**6 + xh**4*(15 - 60*yh**2) + 30*xh**2*yh**2*(-1 + 2*yh**2))*zh)
        gzz=6.831841051919143E-1_r8*(6*ir**2*(xh**6 - 15*xh**4*yh**2 + 15*xh**2*yh**4 - yh**6)*(-1 + 8*zh**2))
    end select
end subroutine

!> calculate the memory used, in bytes
function size_in_mem(basis) result(mem)
    !> basis set
    class(rl_lcao_basis_set), intent(in) :: basis
    !> memory in bytes
    integer :: mem

    integer :: i

    mem=0
    mem=mem+storage_size(basis)
    if ( allocated(basis%species) ) then
        do i=1,size(basis%species)
            mem=mem+storage_size(basis%species(i))
            if ( allocated(basis%species(i)%basis_fn) )      mem=mem+storage_size(basis%species(i)%basis_fn)*size(basis%species(i)%basis_fn)
            if ( allocated(basis%species(i)%basis_l) )       mem=mem+storage_size(basis%species(i)%basis_l)*size(basis%species(i)%basis_l)
            if ( allocated(basis%species(i)%basis_m) )       mem=mem+storage_size(basis%species(i)%basis_m)*size(basis%species(i)%basis_m)
            if ( allocated(basis%species(i)%basis_ang) )     mem=mem+storage_size(basis%species(i)%basis_ang)*size(basis%species(i)%basis_ang)
            if ( allocated(basis%species(i)%basis_group_l) ) mem=mem+storage_size(basis%species(i)%basis_group_l)*size(basis%species(i)%basis_group_l)
            if ( allocated(basis%species(i)%basis_group) )   mem=mem+storage_size(basis%species(i)%basis_group)*size(basis%species(i)%basis_group)
            if ( allocated(basis%species(i)%basis_cutoff) )  mem=mem+storage_size(basis%species(i)%basis_cutoff)*size(basis%species(i)%basis_cutoff)
        enddo
    endif
    if ( allocated(basis%radial) ) then
        do i=1,size(basis%radial)
            mem=mem+storage_size(basis%radial(i))
            if ( allocated(basis%radial(i)%coeff) ) mem=mem+storage_size(basis%radial(i)%coeff)*size(basis%radial(i)%coeff)
            if ( allocated(basis%radial(i)%coeff_deriv) ) mem=mem+storage_size(basis%radial(i)%coeff)*size(basis%radial(i)%coeff_deriv)
        enddo
    endif
    if ( allocated(basis%radial_raw) ) then
        do i=1,size(basis%radial_raw)
            mem=mem+storage_size(basis%radial_raw(i))
            if ( allocated(basis%radial_raw(i)%coeff_deriv) ) mem=mem+storage_size(basis%radial_raw(i)%coeff)*size(basis%radial_raw(i)%coeff_deriv)
        enddo
    endif
    mem=mem/8
end function

!> destroy
subroutine destroy(basis)
    !> basis set
    class(rl_lcao_basis_set), intent(inout) :: basis

    integer :: i
    if ( allocated(basis%species) ) then
        do i=1,size(basis%species)
            if ( allocated(basis%species(i)%basis_fn) )      deallocate(basis%species(i)%basis_fn)
            if ( allocated(basis%species(i)%basis_l) )       deallocate(basis%species(i)%basis_l)
            if ( allocated(basis%species(i)%basis_m) )       deallocate(basis%species(i)%basis_m)
            if ( allocated(basis%species(i)%basis_ang) )     deallocate(basis%species(i)%basis_ang)
            if ( allocated(basis%species(i)%basis_group_l) ) deallocate(basis%species(i)%basis_group_l)
            if ( allocated(basis%species(i)%basis_group) )   deallocate(basis%species(i)%basis_group)
            if ( allocated(basis%species(i)%basis_cutoff) )  deallocate(basis%species(i)%basis_cutoff)
        enddo
        deallocate(basis%species)
    endif
    if ( allocated(basis%radial) ) then
        do i=1,size(basis%radial)
            if ( allocated(basis%radial(i)%coeff) ) deallocate(basis%radial(i)%coeff)
            if ( allocated(basis%radial(i)%coeff_deriv) ) deallocate(basis%radial(i)%coeff_deriv)
        enddo
        deallocate(basis%radial)
    endif
    if ( allocated(basis%radial_raw) ) then
        do i=1,size(basis%radial_raw)
            if ( allocated(basis%radial_raw(i)%coeff) ) deallocate(basis%radial_raw(i)%coeff)
            if ( allocated(basis%radial_raw(i)%coeff_deriv) ) deallocate(basis%radial_raw(i)%coeff_deriv)
        enddo
        deallocate(basis%radial_raw)
    endif
end subroutine

! !> Verbatim copy of the routine from AIMS that evaluates spherical harmonics and gradient, to check.
! subroutine aims_spherical_harm(r,l,m,ylm,dYlmdx,dYlmdy,dYlmdz)
!     real(r8), dimension(3), intent(in) :: r
!     integer, intent(in) :: l,m
!     real(r8), intent(out) :: ylm,dYlmdx,dYlmdy,dYlmdz
!
!
!     real(r8) :: cosph,sinph,costh,sinth
!
!     trig: block
!         real(r8) :: f0,ab,abc
!         ! First use the aims thingy to get the trigonometric guys
!         f0=maxval(abs(r(1:2)))
!         if ( f0 .gt. 0.0_r8 ) then
!             ab = sqrt( r(1)**2 + r(2)**2 )
!             ! cos(phi)
!             cosph = r(1)/ab
!             ! sin(phi)
!             sinph = r(2)/ab
!         else
!             ! we're on the z axis - phi not defined.
!             cosph=1.0_r8
!             sinph=1.0_r8
!             ab = 0.0_r8
!         end if
!
!         ! next, out-of-xy-plane vector theta
!         f0 = max(f0,abs(r(3)) )
!         if ( f0 .gt. 0.0_r8 ) then
!             abc = sqrt( ab**2 + r(3)**2 )
!             ! cos(theta)
!             costh = r(3) / abc
!             ! sin(theta)
!             sinth = ab / abc
!         else
!             ! this piece of code should never trigger, but place us on the z axis anyway
!             costh=1.0_r8
!             sinth=0.0_r8
!         endif
!     end block trig
!
!     ! Get the value of Ylm
!     getval: block
!         real(r8), dimension(3,3) :: m0
!         real(r8), dimension(3) :: v0,v1
!         real(r8), dimension((l+1)**2) :: y,yth,yph
!         real(r8) :: ddth,ddph,f0
!         integer :: i,ll,mm
!
!         !call lo_increment_ylm(sinth,costh,sinph,cosph,0,l,y)
!         call lo_increment_ylm_deriv(sinth,costh,sinph,cosph,0,l,y,yth,yph)
!
!         ylm=-1E10
!         i=0
!         lloop: do ll=0,l
!         do mm=-ll,ll
!             i=i+1
!             if ( ll .eq. l .and. mm .eq. m ) then
!                 ylm=y(i)
!                 ddth=yth(i)
!                 ddph=yph(i)
!                 exit lloop
!             endif
!         enddo
!         enddo lloop
!
!         ! Construct the Cartesian gradient from this?
!         ! this guy should convert from spherical unit vectors to Cartesian.
!         m0(1,:)=[cosph*sinth,costh*cosph,-sinph]
!         m0(2,:)=[sinth*sinph,costh*sinph,cosph]
!         m0(3,:)=[costh,-sinth,0.0_r8]
!         !
!         f0=norm2(r)
!         if ( f0 .gt. 1E-10_r8 ) then
!             v0=[0.0_r8,ddth,ddph]/f0
!         else
!             v0=0.0_r8
!         endif
!         ! To cartesian.
!         v1=matmul(m0,v0)
!         ! And return the gradient
!         dYlmdx=v1(1)
!         dYlmdy=v1(2)
!         dYlmdz=v1(3)
!     end block getval
! end subroutine

! SUBROUTINE lo_increment_ylm_deriv( SINTH, COSTH, SINPH, COSPH, LMIN, LMAX, Y, YTH, YPH)
!     IMPLICIT NONE
!     ! input
!     DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
!     INTEGER            LMIN
!     INTEGER            LMAX
!     ! output
!     ! Y() is recomputed as well, even if it was calculated already earlier.
!     ! This does not cost anything as we need to follow through the entire recursion
!     ! anyhow.
!     ! YTH is dy/d(theta)
!     ! YPH is 1/sin(theta)*[dy/d(phi)]
!     real*8             Y((LMAX+1)*(LMAX+1))
!     real*8             YTH((LMAX+1)*(LMAX+1))
!     real*8             YPH((LMAX+1)*(LMAX+1))
!
!     INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
!     INTEGER            I22L, I24L2
!     DOUBLE PRECISION   D4LL1C, D4LL1S, D2L13, PI
!     DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMP4, TEMP5
!     DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
!     DOUBLE PRECISION   YTHLLR, YTHLL1R, YTHL1L1R, YTHLMR
!     DOUBLE PRECISION   YPHLLR, YPHLL1R, YPHL1L1R, YPHLMR
!     DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!     DOUBLE PRECISION   YTHLLI, YTHLL1I, YTHL1L1I, YTHLMI
!     DOUBLE PRECISION   YPHLLI, YPHLL1I, YPHL1L1I, YPHLMI
!
!     PI = (4.0D+0)*ATAN(1.0D+0)
!
!     if (lmin.le.0) then
!         ! Y(0,0)
!         YLLR = 1.0D+0/SQRT(4.0D+0*PI)
!         YLLI = 0.0D+0
!         Y(1) = YLLR
!         YTH(1) = 0.d0
!         YPH(1) = 0.d0
!         ! continue only if spherical harmonics for (L .GT. 0) are desired
!     end if
!
!     if ( (lmin.le.1).and.(lmax.ge.1)) then
!         !       Y(1,0)
!         Y(3) = SQRT(3.0D+0)*YLLR*COSTH
!         YTH(3) = - SQRT(3.0D+0)*YLLR*SINTH
!         YPH(3) = 0.d0
!         !       Y(1,1) ( = -DCONJG(Y(1,-1)))
!         TEMP1 = SQRT(3.0D+0)*YLLR
!
!         TEMP2 = -TEMP1*SINTH
!         Y(4) = TEMP2*COSPH
!         Y(2) = -TEMP2*SINPH
!         TEMP2 = -TEMP1*COSTH
!         YTH(4) = TEMP2*COSPH
!         YTH(2) = -TEMP2*SINPH
!         YPH(4) = TEMP1*SINPH
!         YPH(2) = TEMP1*COSPH
!     end if
!
!     ! Now calculate remaining ylm's, derivatives in table
!     DO 20 L = max( 2, lmin), LMAX, 1
!         INDEX  = L*L+1
!         INDEX2 = INDEX + 2*L
!         MSIGN  = 1 - 2*MOD(L,2)
!         ! YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!         YL1L1R = Y(INDEX-1)
!         YL1L1I = - MSIGN * Y(INDEX-2*L+1)
!         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))
!         TEMP2 = TEMP1*SINTH
!
!         YLLR = (COSPH*YL1L1R - SINPH*YL1L1I)
!         YLLI = (COSPH*YL1L1I + SINPH*YL1L1R)
!
!         Y(INDEX2) = TEMP2 * YLLR
!         Y(INDEX)  = TEMP2 * MSIGN * YLLI
!
!         ! YTHLL = YTH(L,L) = f(Y(L-1,L-1),YTH(L-1,L-1)) ... Formula 1
!         YTHL1L1R = YTH(INDEX-1)
!         YTHL1L1I = - MSIGN * YTH(INDEX-2*L+1)
!         TEMP3 = TEMP1*COSTH
!
!         YTHLLR = (COSPH*YTHL1L1R - SINPH*YTHL1L1I)
!         YTHLLI = (COSPH*YTHL1L1I + SINPH*YTHL1L1R)
!
!         YTH(INDEX2) = TEMP3 * YLLR + TEMP2 * YTHLLR
!         YTH(INDEX)  = TEMP3 * MSIGN * YLLI + TEMP2 * MSIGN * YTHLLI
!
!         ! YPHLL = YPH(L,L) = f(Y(L-1,L-1),YPH(L-1,L-1)) ... Formula 1
!         YPHL1L1R = YPH(INDEX-1)
!         YPHL1L1I = - MSIGN * YPH(INDEX-2*L+1)
!
!         YLLR = (-SINPH*YL1L1R - COSPH*YL1L1I)
!         YLLI = (-SINPH*YL1L1I + COSPH*YL1L1R)
!
!         YPHLLR = (COSPH*YPHL1L1R - SINPH*YPHL1L1I)
!         YPHLLI = (COSPH*YPHL1L1I + SINPH*YPHL1L1R)
!
!         YPH(INDEX2) = TEMP1 * YLLR + TEMP2 * YPHLLR
!         YPH(INDEX)  = TEMP1 * MSIGN * YLLI + TEMP2 * MSIGN * YPHLLI
!         !
!         INDEX2 = INDEX2 - 1
!         INDEX  = INDEX  + 1
!         !
!         !        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!         !               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!         !
!         TEMP1 = SQRT(DBLE(2*L+1))
!         TEMP2 = TEMP1*COSTH
!         YLL1R = TEMP2*YL1L1R
!         YLL1I = TEMP2*YL1L1I
!         Y(INDEX2) = YLL1R
!         Y(INDEX)  = - MSIGN * YLL1I
!         !
!         !        YTHLL1 = YTH(L,L-1) = f(Y(L-1,L-1),YTH(L-1,L-1)) ... Formula 2
!         !               (the coefficient for Y(L-2,L-1),YTH(L-2,L-1) in Formula 2 is zero)
!         !
!         TEMP3 = -TEMP1*SINTH
!         YLL1R = TEMP3*YL1L1R
!         YLL1I = TEMP3*YL1L1I
!         YTHLL1R = TEMP2*YTHL1L1R
!         YTHLL1I = TEMP2*YTHL1L1I
!         YTH(INDEX2) = YLL1R + YTHLL1R
!         YTH(INDEX)  = - MSIGN * ( YLL1I + YTHLL1I )
!
!         !
!         !        YPHLL1 = YPH(L,L-1) = f(YPH(L-1,L-1)) ... Formula 2
!         !               (the coefficient for YPH(L-2,L-1) in Formula 2 is zero)
!         !
!         !         TEMP2 = TEMP2*SINTH
!         YPHLL1R = TEMP2*YPHL1L1R
!         YPHLL1I = TEMP2*YPHL1L1I
!         YPH(INDEX2) = YPHLL1R
!         YPH(INDEX)  = - MSIGN * YPHLL1I
!
!         INDEX2 = INDEX2 - 1
!         INDEX  = INDEX  + 1
!         !
!         I4L2 = INDEX - 4*L + 2
!         I2L  = INDEX - 2*L
!         I24L2 = INDEX2 - 4*L + 2
!         I22L  = INDEX2 - 2*L
!         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
!         D4LL1S = SINTH*SQRT(DBLE(4*L*L-1))
!         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!         !
!         DO 10 M = L - 2, 0, -1
!             !
!             !        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!             !
!             TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
!
!             TEMP2 = D4LL1C*TEMP1
!             TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
!             YLMR = TEMP2*Y(I22L) + TEMP3*Y(I24L2)
!             YLMI = TEMP2*Y(I2L) + TEMP3*Y(I4L2)
!             Y(INDEX2) = YLMR
!             Y(INDEX)  = YLMI
!             !
!             TEMP4 = - D4LL1S*TEMP1
!             YTHLMR = TEMP4*Y(I22L) + TEMP2*YTH(I22L) + TEMP3*YTH(I24L2)
!             YTHLMI = TEMP4*Y(I2L) + TEMP2*YTH(I2L) + TEMP3*YTH(I4L2)
!             YTH(INDEX2) = YTHLMR
!             YTH(INDEX)  = YTHLMI
!             !
!             TEMP2 = D4LL1C*TEMP1
!             TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
!             YPHLMR = TEMP2*YPH(I22L) + TEMP3*YPH(I24L2)
!             YPHLMI = TEMP2*YPH(I2L) + TEMP3*YPH(I4L2)
!             YPH(INDEX2) = YPHLMR
!             YPH(INDEX)  = YPHLMI
!             !
!             INDEX2 = INDEX2 - 1
!             INDEX  = INDEX  + 1
!             I24L2   = I24L2   - 1
!             I22L    = I22L    - 1
!             I4L2   = I4L2   + 1
!             I2L    = I2L    + 1
!         10    CONTINUE
!     20 CONTINUE
! end subroutine

! subroutine lo_increment_ylm( sinth, costh, sinph, cosph, lmin, lmax, y)
!     real(r8) :: costh, sinth, cosph, sinph
!     integer, intent(in) :: lmin,lmax
!     real(r8) :: y( (lmax+1)*(lmax+1) )
!
!     !  VB, 10/08/05:
!     !  modified version of original subroutine; spherical angles now passed directly as
!     !  argument, and it is possible to simply increment a table of already given ylm
!     !
!     !  VB, 11/23/04:
!     !  Subroutine ylm obtained from Ricardo Gomez-Abral, the original source is not
!     !  known to me.
!     !  Given the absence of a formal copyright notice, I assume that the code is in the
!     !  public domain; however, I do not understand the signature "Institut fuer theoretische
!     !  Chemie - TU Vienna".
!     !
!     !  VB, 07.02.05:
!     !  created this modified version which obtains the real Ylm only, i.e. does not
!     !  bake them together to a complex version as does the original ylm.f
!
!     !     ..................................................................
!     ! 1.     PROGRAM UNIT 'YLM'
!     !           Calculates spherical harmonics
!     !           FORTRAN 77 SUBROUTINE
!     !
!     ! 2.     PURPOSE
!     !           The spherical harmonics (Condon and Shortley convention)
!     !             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!     !           for vector V (given in Cartesian coordinates)
!     !           are calculated. In the Condon Shortley convention the
!     !           spherical harmonics are defined as
!     !                             +------+
!     !                        m    |   1     m              im(Phi)
!     !           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!     !                            \| 2(Pi)   l
!     !                  m
!     !           where P (cos(Theta)) is the normalized Associated Legendre
!     !                  l
!     !           function. Thus,
!     !                                          m      *
!     !                            Y(l,-m) = (-1) Y(l,m)
!     !
!     !
!     ! 3.     USAGE
!     !           DOUBLE PRECISION V(3), Y(5*5)
!     !           V(1) = ...
!     !           V(2) = ...
!     !           V(3) = ...
!     !           CALL YLM(V,4,Y)
!     !
!     !        ARGUMENT-DESCRIPTION
!     !           SINPH, COSPH, SINTH, COSTH - DOUBLE PRECISION variables        (input)
!     !                    the
!     !                    angles Theta and Phi necessary for the calculation
!     !                    of the spherical harmonics.
!     !           LMIN   - INTEGER value                               (input)
!     !                    lower bound, below which spherical harmonics are already given.
!     !                    if LMIN > 0, expect that Y(0 ... LMIN-1) are already tabulated.
!     !           LMAX   - INTEGER value                               (input)
!     !                    upper bound of L for which spherical harmonics
!     !                    will be calculated
!     !                    constraint:
!     !                       LMAX .GE. 0 (not checked)
!     !           Y      - REAL*8 array, dimension (LMAX+1)**2    (output)
!     !                    contains the calculated spherical harmonics
!     !                    Y(1)                   for L .EQ. 0 (M = 0)
!     !                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!     !                    ...
!     !                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!     !                                           for L .EQ. LMAX
!     !                                              (M = -L,...,L)
!     !                    constraint:
!     !                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!     !
!     !        USED SUBROUTINES (DIRECTLY CALLED)
!     !           none
!     !
!     !        INDIRECTLY CALLED SUBROUTINES
!     !           none
!     !
!     !        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!     !           none
!     !
!     !        INPUT/OUTPUT (READ/WRITE)
!     !           none
!     !
!     !        MACHINENDEPENDENT PROGRAMPARTS
!     !           Type COMPLEX*16 is used which does not conform to the
!     !           FORTRAN 77 standard.
!     !           Also the non-standard type conversion function DCMPLX()
!     !           is used which combines two double precision values into
!     !           one double complex value.
!     !
!     ! 4.     REMARKS
!     !           none
!     !
!     ! 5.     METHOD
!     !           The basic algorithm used to calculate the spherical
!     !           harmonics for vector V is as follows:
!     !
!     !           Y(0,0)
!     !           Y(1,0)
!     !           Y(1,1)
!     !           Y(1,-1) = -Y(1,1)
!     !           DO L = 2, LMAX
!     !              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!     !              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!     !              DO M = L-2, 0, -1
!     !                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!     !                 Y(L,-M)= (-1)**M*Y(L,M)
!     !              ENDDO
!     !           ENDDO
!     !
!     !           In the following the necessary recursion formulas and
!     !           starting values are given:
!     !
!     !        Start:
!     !                        +------+
!     !                        |   1
!     !           Y(0,0) =  -+ | -----
!     !                       \| 4(Pi)
!     !
!     !                                   +------+
!     !                                   |   3
!     !           Y(1,0) =  cos(Theta) -+ | -----
!     !                                  \| 4(Pi)
!     !
!     !                                     +------+
!     !                                     |   3    i(Phi)
!     !           Y(1,1) =  - sin(Theta) -+ | ----- e
!     !                                    \| 8(Pi)
!     !
!     !        Formula 1:
!     !
!     !           Y(l,l) =
!     !                           +--------+
!     !                           | (2l+1)   i(Phi)
!     !            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!     !                          \|   2l
!     !
!     !        Formula 2:
!     !                                  +---------------+
!     !                                  |  (2l-1)(2l+1)
!     !           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!     !                                 \|   (l-m)(l+m)
!     !
!     !                                    +--------------------+
!     !                                    |(l-1+m)(l-1-m)(2l+1)
!     !                              -  -+ |-------------------- Y(l-2,m)
!     !                                   \|  (2l-3)(l-m)(l+m)
!     !
!     !        Formula 3: (not used in the algorithm because of the division
!     !                    by sin(Theta) which may be zero)
!     !
!     !                                    +--------------+
!     !                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!     !           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!     !                      sin(Theta)   \| (l+m+1)(l-m)
!     !
!     !                                    +--------------+
!     !                                    |(l-m-1)(l+m+2)  -2i(Phi)
!     !                              -  -+ |-------------- e        Y(l,m+2)
!     !                                   \| (l-m)(l+m+1)
!     !
!     !
!     !        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     !     ..................................................................
!     !
!     !
!     integer :: i2l, i4l2, index, index2, l, m, msign
!     !vb
!     integer :: i22l, i24l2
!     !vb end
!     real(r8) :: d4ll1c, d2l13, pi
!     real(r8) :: temp1, temp2, temp3
!     real(r8) :: yllr, yll1r, yl1l1r, ylmr
!     real(r8) :: ylli, yll1i, yl1l1i, ylmi
!
!     pi = (4.0d+0)*atan(1.0d+0)
!
!     if (lmin.le.0) then
!         ! Y(0,0)
!          YLLR = 1.0D+0/SQRT(4.0D+0*PI)
!          YLLI = 0.0D+0
!          Y(1) = YLLR
!         ! continue only if spherical harmonics for (L .GT. 0) are desired
!     end if
!
!     if ( (lmin.le.1).and.(lmax.ge.1)) then
!         ! Y(1,0)
!         Y(3) = SQRT(3.0D+0)*YLLR*COSTH
!         ! Y(1,1) ( = -DCONJG(Y(1,-1)))
!         TEMP1 = -SQRT(3.0D+0)*YLLR*SINTH
!         Y(4) = TEMP1*COSPH
!         Y(2) = -TEMP1*SINPH
!     end if
!
!     ! now calculate remaining ylm's in table
!     do l = max( 2, lmin), lmax
!         index  = l*l+1
!         index2 = index + 2*l
!         msign  = 1 - 2*mod(l,2)
!         ! yll = y(l,l) = f(y(l-1,l-1)) ... formula 1
!         yl1l1r = y(index-1)
!         yl1l1i = - msign * y(index-2*l+1)
!         temp1 = -sqrt(dble(2*l+1)/dble(2*l))*sinth
!         yllr = temp1*(cosph*yl1l1r - sinph*yl1l1i)
!         ylli = temp1*(cosph*yl1l1i + sinph*yl1l1r)
!         y(index2) = yllr
!         y(index)  = msign * ylli
!         index2 = index2 - 1
!         index  = index  + 1
!         ! yll1 = y(l,l-1) = f(y(l-1,l-1)) ... formula 2
!         ! (the coefficient for y(l-2,l-1) in formula 2 is zero)
!         temp2 = sqrt(dble(2*l+1))*costh
!         yll1r = temp2*yl1l1r
!         yll1i = temp2*yl1l1i
!         y(index2) = yll1r
!         y(index)  = - msign * yll1i
!         index2 = index2 - 1
!         index  = index  + 1
!
!         i4l2 = index - 4*l + 2
!         i2l  = index - 2*l
!         i24l2 = index2 - 4*l + 2
!         i22l  = index2 - 2*l
!         d4ll1c = costh*sqrt(dble(4*l*l-1))
!         d2l13  = -sqrt(dble(2*l+1)/dble(2*l-3))
!         do m = l - 2, 0, -1
!             ! ylm = y(l,m) = f(y(l-2,m),y(l-1,m)) ... formula 2
!             temp1 = 1.0d+0/sqrt(dble((l+m)*(l-m)))
!             temp2 = d4ll1c*temp1
!             temp3 = d2l13*sqrt(dble((l+m-1)*(l-m-1)))*temp1
!             ylmr = temp2*y(i22l) + temp3*y(i24l2)
!             ylmi = temp2*y(i2l) + temp3*y(i4l2)
!             y(index2) = ylmr
!             y(index)  = ylmi
!             !
!             index2 = index2 - 1
!             index  = index  + 1
!             i24l2   = i24l2   - 1
!             i22l    = i22l    - 1
!             i4l2   = i4l2   + 1
!             i2l    = i2l    + 1
!         enddo
!     enddo
! end subroutine

! !> Central difference rules up to pretty high order
! subroutine lo_centraldifference(n,x0,h,sc)
!     !> order of stencil
!     integer, intent(in) :: n
!     !> point to differentiate around
!     real(r8), intent(in) :: x0
!     !> step
!     real(r8), intent(in) :: h
!     !> abscissas and weights. Should be dimension (2*n+1,4), where the columns are abscissa, and weights for 1st,2nd and 3rd derivatives
!     real(r8), dimension(:,:), intent(out) :: sc
!
!     sc=0.0_r8
!     select case(n)
!     case(1)
!         sc(1,:)=[-1.0_r8,-0.50000000000000_r8, 1.0000000000000_r8,0.0_r8];
!         sc(2,:)=[ 0.0_r8, 0.0_r8,              -2.0000000000000_r8,0.0_r8];
!         sc(3,:)=[ 1.0_r8, 0.50000000000000_r8, 1.0000000000000_r8,0.0_r8];
!     case(2)
!         sc(1,:)=[-2.0_r8, 8.333333333333333e-2_r8, -8.333333333333333e-2_r8, -5.000000000000000e-1_r8];
!         sc(2,:)=[-1.0_r8, -6.666666666666667e-1_r8, 1.333333333333333_r8, 1.000000000000000_r8];
!         sc(3,:)=[0.0_r8, 0.0_r8, -2.500000000000000_r8, 0.0_r8];
!         sc(4,:)=[1.0_r8, 6.666666666666667e-1_r8, 1.333333333333333_r8, -1.000000000000000_r8];
!         sc(5,:)=[2.0_r8, -8.333333333333333e-2_r8, -8.333333333333333e-2_r8, 5.000000000000000e-1_r8];
!     case(3)
!         sc(1,:)=[-3.0_r8, -1.666666666666667e-2_r8, 1.111111111111111e-2_r8, 1.250000000000000e-1_r8];
!         sc(2,:)=[-2.0_r8, 1.500000000000000e-1_r8, -1.500000000000000e-1_r8, -1.000000000000000_r8];
!         sc(3,:)=[-1.0_r8, -7.500000000000000e-1_r8, 1.500000000000000_r8, 1.625000000000000_r8];
!         sc(4,:)=[0.0_r8, 0.0_r8, -2.722222222222222_r8, 0.0_r8];
!         sc(5,:)=[1.0_r8, 7.500000000000000e-1_r8, 1.500000000000000_r8, -1.625000000000000_r8];
!         sc(6,:)=[2.0_r8, -1.500000000000000e-1_r8, -1.500000000000000e-1_r8, 1.000000000000000_r8];
!         sc(7,:)=[3.0_r8, 1.666666666666667e-2_r8, 1.111111111111111e-2_r8, -1.250000000000000e-1_r8];
!     case(4)
!         sc(1,:)=[-4.0_r8, 3.571428571428571e-3_r8, -1.785714285714286e-3_r8, -2.916666666666667e-2_r8];
!         sc(2,:)=[-3.0_r8, -3.809523809523810e-2_r8, 2.539682539682540e-2_r8, 3.000000000000000e-1_r8];
!         sc(3,:)=[-2.0_r8, 2.000000000000000e-1_r8, -2.000000000000000e-1_r8, -1.408333333333333_r8];
!         sc(4,:)=[-1.0_r8, -8.000000000000000e-1_r8, 1.600000000000000_r8, 2.033333333333333_r8];
!         sc(5,:)=[0.0_r8, 0.0_r8, -2.847222222222222_r8, 0.0_r8];
!         sc(6,:)=[1.0_r8, 8.000000000000000e-1_r8, 1.600000000000000_r8, -2.033333333333333_r8];
!         sc(7,:)=[2.0_r8, -2.000000000000000e-1_r8, -2.000000000000000e-1_r8, 1.408333333333333_r8];
!         sc(8,:)=[3.0_r8, 3.809523809523810e-2_r8, 2.539682539682540e-2_r8, -3.000000000000000e-1_r8];
!         sc(9,:)=[4.0_r8, -3.571428571428571e-3_r8, -1.785714285714286e-3_r8, 2.916666666666667e-2_r8];
!     case(5)
!         sc(1,:)=[-5.0_r8, -7.936507936507937e-4_r8, 3.174603174603175e-4_r8, 6.779100529100529e-3_r8];
!         sc(2,:)=[-4.0_r8, 9.920634920634921e-3_r8, -4.960317460317460e-3_r8, -8.339947089947090e-2_r8];
!         sc(3,:)=[-3.0_r8, -5.952380952380952e-2_r8, 3.968253968253968e-2_r8, 4.830357142857143e-1_r8];
!         sc(4,:)=[-2.0_r8, 2.380952380952381e-1_r8, -2.380952380952381e-1_r8, -1.733730158730159_r8];
!         sc(5,:)=[-1.0_r8, -8.333333333333333e-1_r8, 1.666666666666667_r8, 2.318055555555556_r8];
!         sc(6,:)=[0.0_r8, 0.0_r8, -2.927222222222222_r8, 0.0_r8];
!         sc(7,:)=[1.0_r8, 8.333333333333333e-1_r8, 1.666666666666667_r8, -2.318055555555556_r8];
!         sc(8,:)=[2.0_r8, -2.380952380952381e-1_r8, -2.380952380952381e-1_r8, 1.733730158730159_r8];
!         sc(9,:)=[3.0_r8, 5.952380952380952e-2_r8, 3.968253968253968e-2_r8, -4.830357142857143e-1_r8];
!         sc(10,:)=[4.0_r8, -9.920634920634921e-3_r8, -4.960317460317460e-3_r8, 8.339947089947090e-2_r8];
!         sc(11,:)=[5.0_r8, 7.936507936507937e-4_r8, 3.174603174603175e-4_r8, -6.779100529100529e-3_r8];
!     case(6)
!         sc(1,:)=[-6.0_r8, 1.803751803751804e-4_r8, -6.012506012506013e-5_r8, -1.583994708994709e-3_r8];
!         sc(2,:)=[-5.0_r8, -2.597402597402597e-3_r8, 1.038961038961039e-3_r8, 2.261904761904762e-2_r8];
!         sc(3,:)=[-4.0_r8, 1.785714285714286e-2_r8, -8.928571428571429e-3_r8, -1.530952380952381e-1_r8];
!         sc(4,:)=[-3.0_r8, -7.936507936507937e-2_r8, 5.291005291005291e-2_r8, 6.572751322751323e-1_r8];
!         sc(5,:)=[-2.0_r8, 2.678571428571429e-1_r8, -2.678571428571429e-1_r8, -1.995089285714286_r8];
!         sc(6,:)=[-1.0_r8, -8.571428571428571e-1_r8, 1.714285714285714_r8, 2.527142857142857_r8];
!         sc(7,:)=[0.0_r8, 0.0_r8, -2.982777777777778_r8, 0.0_r8];
!         sc(8,:)=[1.0_r8, 8.571428571428571e-1_r8, 1.714285714285714_r8, -2.527142857142857_r8];
!         sc(9,:)=[2.0_r8, -2.678571428571429e-1_r8, -2.678571428571429e-1_r8, 1.995089285714286_r8];
!         sc(10,:)=[3.0_r8, 7.936507936507937e-2_r8, 5.291005291005291e-2_r8, -6.572751322751323e-1_r8];
!         sc(11,:)=[4.0_r8, -1.785714285714286e-2_r8, -8.928571428571429e-3_r8, 1.530952380952381e-1_r8];
!         sc(12,:)=[5.0_r8, 2.597402597402597e-3_r8, 1.038961038961039e-3_r8, -2.261904761904762e-2_r8];
!         sc(13,:)=[6.0_r8, -1.803751803751804e-4_r8, -6.012506012506013e-5_r8, 1.583994708994709e-3_r8];
!     end select
!     ! Now adjust it to the interval we care about
!     sc(:,1)=sc(:,1)*h+x0
!     sc(:,2)=sc(:,2)/h
!     sc(:,3)=sc(:,3)/h**2
!     sc(:,4)=sc(:,4)/h**3
! end subroutine

end module
