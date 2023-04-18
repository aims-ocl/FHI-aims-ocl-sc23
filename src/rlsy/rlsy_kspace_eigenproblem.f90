module rlsy_kspace_eigenproblem
!!
!! Handle that deals with the Kohn-Sham eigenproblem. Stores the relevant parts of the
!! solution. There are two ways this can be rearranged, that will make the code flow
!! quite differently, so it's worth trying to understand them properly:
!!
!! Single-proc mode
!!
!! ┌────┐    ┌───────────────────────┐       ┌─────────────────────────────────────────────┐
!! │    ├────┤  k-point 1 on rank 0  ├───┬───┤KS eigenvalues (n_basis x n_spin)            │
!! │    │    └───────────────────────┘   │   ├─────────────────────────────────────────────┤
!! │MPI │    ┌───────────────────────┐   ├───┤KS eigenvectors (n_basis x n_basis x n_spin) │
!! │Rank├────┤  k-point 2 on rank 0  │   │   ├─────────────────────────────────────────────┤
!! │ 0  │    └───────────────────────┘   └───┤Density matrix (n_basis x n_basis x n_spin)  │
!! │    │    ┌───────────────────────┐       └─────────────────────────────────────────────┘
!! │    ├────┤  k-point 3 on rank 0  │
!! └────┘    └───────────────────────┘
!! ┌────┐    ┌───────────────────────┐
!! │    ├────┤  k-point 1 on rank 1  │
!! │    │    └───────────────────────┘
!! │MPI │    ┌───────────────────────┐
!! │Rank├────┤  k-point 2 on rank 1  │
!! │ 1  │    └───────────────────────┘
!! │    │    ┌───────────────────────┐
!! │    ├────┤  k-point 3 on rank 1  │
!! └────┘    └───────────────────────┘
!!
!! Multi-proc mode
!!
!! ┌───────┐  ┌──────────────────────────┐  ┌──────────────┐
!! │       ├──┤MPI Rank 0 local, 0 global├──┤              │
!! │k-point│  └──────────────────────────┘  │Density matrix│
!! │   1   │  ┌──────────────────────────┐  │              │
!! │ spin  ├──┤MPI Rank 1 local, 1 global├──┤ distributed  │
!! │channel│  └──────────────────────────┘  │ across ranks │
!! │   1   │  ┌──────────────────────────┐  │              │
!! │       ├──┤MPI Rank 2 local, 2 global├──┤              │
!! └───────┘  └──────────────────────────┘  └──────────────┘
!! ┌───────┐  ┌──────────────────────────┐  ┌──────────────┐
!! │       ├──┤MPI Rank 0 local, 3 global├──┤              │
!! │k-point│  └──────────────────────────┘  │Density matrix│
!! │   1   │  ┌──────────────────────────┐  │              │
!! │ spin  ├──┤MPI Rank 1 local, 4 global├──┤ distributed  │
!! │channel│  └──────────────────────────┘  │ across ranks │
!! │   2   │  ┌──────────────────────────┐  │              │
!! │       ├──┤MPI Rank 2 local, 5 global├──┤              │
!! └───────┘  └──────────────────────────┘  └──────────────┘
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_exitcode_param,rl_iou,&
    rl_exitcode_memory,rl_exitcode_symmetry,rl_sqtol,rl_exitcode_mpi
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: rl_chop,tochar,rl_clean_fractional_coordinates,rl_sqnorm
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_sorting, only: rl_return_unique
use rlsy_verletlist, only: rl_verletbox
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_scalapack_helper, only: rl_blacs_helper
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_kpointmesh, only: rl_kpoint_mesh
use elsi !, only: elsi_handle,elsi_init

implicit none
private
public :: rl_kspace_eigenproblem
public :: rl_kspace_eigenproblem_singleproc
public :: rl_kspace_eigenproblem_multiproc
public :: rl_kspace_eigenproblem_prob
public :: rl_kspace_eigenproblem_prob_real
public :: rl_kspace_eigenproblem_prob_complex
public :: rl_kspace_eigenproblem_kpoint_real
public :: rl_kspace_eigenproblem_kpoint_complex
public :: rl_setup_kspace_eigenproblem

! Some parameters for ELSI
real(r8), parameter :: zero_threshold=1E-15_r8
integer, parameter :: solver_ELPA=1

!> common information about one k-point
type :: rl_kspace_eigenproblem_kpoint
    !> irreducible index to this k-point
    integer :: irreducible_index=-rl_hugeint
    !> KS eigenvalues (n_basis,n_spin)
    real(r8), dimension(:,:), allocatable :: eigenvalue
    !> Occupation number (n_state,n_spin)
    real(r8), dimension(:,:), allocatable :: occupation
end type
type, extends(rl_kspace_eigenproblem_kpoint) :: rl_kspace_eigenproblem_kpoint_real
    !> buffer for eigenvectors (n_basis,n_basis,n_spin)
    real(r8), dimension(:,:,:), allocatable :: eigenvector
    !> buffer for density matrix (n_basis,n_basis,n_spin)
    real(r8), dimension(:,:,:), allocatable :: densitymatrix
end type
type, extends(rl_kspace_eigenproblem_kpoint) :: rl_kspace_eigenproblem_kpoint_complex
    !> buffer for eigenvectors (n_basis,n_basis,n_spin)
    complex(r8), dimension(:,:,:), allocatable :: eigenvector
    !> buffer for density matrix (n_basis,n_basis,n_spin)
    complex(r8), dimension(:,:,:), allocatable :: densitymatrix
end type

!> Common things for a single distributed eigenproblem.
type, abstract :: rl_kspace_eigenproblem_prob
    !> local number of rows
    integer :: n_row_local=-rl_hugeint
    !> local number of columns
    integer :: n_col_local=-rl_hugeint
    !> blocksize
    integer :: blocksize=-rl_hugeint
    !> buffer for Kohn-Sham eigenvalues
    real(r8), dimension(:), allocatable :: eigenvalue
end type
type, extends(rl_kspace_eigenproblem_prob) :: rl_kspace_eigenproblem_prob_real
    !> buffer for local Hamiltonian @TODO retire, probably
    real(r8), dimension(:,:), allocatable :: hamiltonian
    !> buffer for local overlap @TODO retire, probably
    real(r8), dimension(:,:), allocatable :: overlap
    !> buffer for eigenvectors
    real(r8), dimension(:,:), allocatable :: eigenvector
    !> buffer for density matrix
    real(r8), dimension(:,:), allocatable :: densitymatrix
end type
type, extends(rl_kspace_eigenproblem_prob) :: rl_kspace_eigenproblem_prob_complex
    !> buffer for local Hamiltonian @TODO retire, probably
    complex(r8), dimension(:,:), allocatable :: hamiltonian
    !> buffer for local overlap @TODO retire, probably
    complex(r8), dimension(:,:), allocatable :: overlap
    !> buffer for eigenvectors
    complex(r8), dimension(:,:), allocatable :: eigenvector
    !> buffer for density matrix
    complex(r8), dimension(:,:), allocatable :: densitymatrix
end type

!> handle that sorts out how to solve the generalized eigenvalue problem
type, abstract :: rl_kspace_eigenproblem
    !> how many states are we interested in?
    integer :: n_state=-rl_hugeint
    !> number of basis functions per cell
    integer :: n_basis=-rl_hugeint
    !> number of electrons
    real(r8) :: n_electron=-rl_huge
    !> number of spin channels
    integer :: n_spin=-rl_hugeint
    !> chemical potential?
    real(r8) :: chemical_potential=-rl_huge
    !> electron entropy
    real(r8) :: entropy=-rl_huge
    !> ELSI handle
    type(elsi_handle) :: eh
end type

!> this variant handles when all eigenproblems are solved simultaneously by many ranks.
type, extends(rl_kspace_eigenproblem) :: rl_kspace_eigenproblem_multiproc
    !> how many groups of mpi ranks are there in total?
    integer :: n_group_global=-rl_hugeint
    !> which group of ranks does this rank belong to?
    integer :: group_index=-rl_hugeint
    !> mpi communicator for this group.
    type(rl_mpi_helper) :: ml
    !> BLACS context for this group
    type(rl_blacs_helper) :: bl
    !> irreducible index of the relevant k-point
    integer :: index_irreducible_kpoint=-rl_hugeint
    !> spin index
    integer :: index_spin=-rl_hugeint
    !> space for local buffers
    class(rl_kspace_eigenproblem_prob), allocatable :: buf
end type

!> this variant is for the case when there are more k-points than ranks, or at least comparable numbers
type, extends(rl_kspace_eigenproblem) :: rl_kspace_eigenproblem_singleproc
    !> how many k-points does this rank handle?
    integer :: n_kpoint
    !> handle for that k-point
    class(rl_kspace_eigenproblem_kpoint), dimension(:), allocatable :: kpoint
end type

contains

!> Set up everything needed to solve the KS equations with ELSI.
subroutine rl_setup_kspace_eigenproblem(KS,p,basis,kmesh,n_states,n_spin,n_electrons,mw,mem,verbosity,&
                elsi_solver,elsi_out_level,occupation_type,solver_method,elsi_elpa_n_single,&
                occupation_width,occupation_acc,basis_threshold)
    !> container for eigenvalues and eigenvectors and such
    class(rl_kspace_eigenproblem), allocatable, intent(out) :: KS
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> k-point mesh
    type(rl_kpoint_mesh), intent(in) :: kmesh
    !> how many states to solve for
    integer, intent(in) :: n_states
    !> how many spin channels
    integer, intent(in) :: n_spin
    !> how many electrons
    real(r8), intent(in) :: n_electrons
    !> MPI helper (global)
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> these are ELSI/ELPA settings
    integer, intent(in) :: elsi_solver
    integer, intent(in) :: elsi_out_level
    integer, intent(in) :: occupation_type
    integer, intent(in) :: solver_method
    integer, intent(in) :: elsi_elpa_n_single
    real(r8), intent(in) :: occupation_width
    real(r8), intent(in) :: occupation_acc
    real(r8), intent(in) :: basis_threshold

    real(r8) :: timer,t0,t1
    logical :: complex_eigenvectors
    !basic: block
        integer, dimension(:), allocatable :: di
        integer :: i,j,l,nprob
        integer :: parallel_mode
    !realorcomplex: block
        real(r8), dimension(3) :: v0
        !integer :: i,l
    !sng: block
        type(rl_mpi_helper) :: ml
        type(rl_blacs_helper) :: bl
        !integer :: i,l
    !mlt: block
        integer, dimension(:,:), allocatable :: probind
        !integer, dimension(:), allocatable :: di,dj
        integer, dimension(:), allocatable :: dj
        !integer :: nprob,i,j,k,l
        integer :: k

    ! store some basic information
    !basic: block

        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
            t0=timer
            t1=timer
            write(rl_iou,*) ''
            write(rl_iou,*) 'SETTING UP EIGENVALUE SOLVER'
        endif

        ! Number of KS problems I have to solve, eventually
        nprob=kmesh%n_irr_kpoint*n_spin

        ! Some dummy arrays to help with distribution
        call mem%allocate(di,nprob,persistent=.false.,scalable=.false.)
        di=0

        ! Count eigenvalue problems per rank
        do i=1,nprob
        do j=1,mw%n
            if ( mod(j,nprob) .eq. i-1) di(i)=di(i)+1
        enddo
        enddo

        ! This should be enough to decide how we are going to solve things.
        ! There is likely no point in using parallel solver unless there are
        ! at least two ranks per problem, otherwise we will just wait a lot,
        ! I suppose.
        if ( minval(di) .lt. 2 ) then
            ! One rank per problem.
            parallel_mode=1
            allocate(rl_kspace_eigenproblem_singleproc::KS)
        else
            ! Many ranks per problem
            parallel_mode=0
            allocate(rl_kspace_eigenproblem_multiproc::KS)
        endif

        ! This should be all of it, but there are some corner cases we have to think about.
        if ( parallel_mode == 1 ) then
            i=mw%n/nprob ! smallest number of ranks per problem.
            if ( KS%n_basis**2 .lt. i ) then
                call rl_stop_gracefully([&
                'There are more ranks per eigenvalue problem than there are ',&
                'components in the Hamiltonian. This is highly inefficient  ',&
                'and most likely a symptom of strange input/runtime choices.']&
                ,rl_exitcode_mpi,mw%comm)
            endif
        endif

        !@TODO Have a hard switch, if the size of the Hamiltonian is less
        ! than some number, always do it single-proc. But for that I need
        ! to know what the crossover is. Maybe a waste of time.

        ! We can set the common things that are always the same:
        KS%n_state = n_states
        KS%n_spin = n_spin
        KS%n_electron = n_electrons
        KS%n_basis = basis%n_basis

        ! A little cleanup
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    !end block basic

    ! Check if we want real or complex eigenvectors
    !realorcomplex: block

        l=0
        kpl: do i=1,kmesh%n_irr_kpoint
            ! Check if k-point real or complex?
            v0=kmesh%ip(i)%r
            v0=matmul(p%inv_reciprocal_latticevectors,v0)
            v0=rl_clean_fractional_coordinates(v0+v0+0.5_r8)-0.5_r8
            v0=rl_chop(v0,1E-11_r8)
            v0=rl_clean_fractional_coordinates(v0+v0+0.5_r8)-0.5_r8
            v0=rl_chop(v0,1E-11_r8)
            if ( rl_sqnorm(v0) .lt. rl_sqtol ) then
                ! k-vector equal to negative itself, real eigenvectors!
                l=l+1
            else
                ! k-vector not equal to negative itsels, complex eigenvectors!
                exit
            endif
        enddo kpl

        if ( l .eq. kmesh%n_irr_kpoint ) then
            ! All of them real!
            complex_eigenvectors=.false.
        else
            ! Some not real. Make all complex.
            complex_eigenvectors=.true.
        endif

        !@TODO Could be smart to choose this on a k-by-k basis, then assign
        ! slightly more ranks to the complex points. Should save a little time.
    !end block realorcomplex

    ! So I made it an abstract class, the handle for the eigenproblem. I think that makes
    ! somewhat sense. It means you can not accidentally try to use the wrong things in the
    ! wrong context, since it won't even be accessible. But the flow of the code, except at
    ! the stage of solving the KS equations should be identical.
    select type(KS)
    type is(rl_kspace_eigenproblem_singleproc)
    ! Here I set it up for processing in the ELSI single-proc mode. I don't divide
    ! between spin channels and stuff. Why does ELSI not support that? I don't know.
    !sng: block

        ! Count k-points per rank
        KS%n_kpoint=0
        do i=1,kmesh%n_irr_kpoint
            if ( mod(i,mw%n) .eq. mw%r ) KS%n_kpoint=KS%n_kpoint+1
        enddo
        ! Make space for k-points
        if ( complex_eigenvectors ) then
            allocate(rl_kspace_eigenproblem_kpoint_complex::KS%kpoint(KS%n_kpoint))
        else
            allocate(rl_kspace_eigenproblem_kpoint_real::KS%kpoint(KS%n_kpoint))
        endif

        ! Store at least some information about the k-point. More will likely come.
        l=0
        do i=1,kmesh%n_irr_kpoint
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            l=l+1
            KS%kpoint(l)%irreducible_index=i
        enddo

        ! Initialize ELSI!
        call elsi_init(KS%eh,&
            solver=solver_ELPA,&
            parallel_mode=0,&
            matrix_format=0,&
            n_basis=KS%n_basis,&
            n_electron=KS%n_electron,&
            n_state=KS%n_state)

        ! In single-proc mode, we use some kind of fake matrix storage thing
        ! split the communicator into single-rank communicators
        call mw%split(ml,mw%r)
        ! create single-communicator blacs context
        call bl%init(ml)
        ! then set blacs storage, with the blocksize equal to n_basis
        call elsi_set_blacs(KS%eh,bl%icontxt,KS%n_basis)
        !@TODO Figure out if it's safe to destroy blacs context and MPI communicator now. Not sure.

        ! Now make space for buffers? Yes no maybe. In single-proc mode, I work under the reasonable
        ! assumption that memory is not an issue, which should hold pretty much all the time. It would
        ! be some rather odd corner-case of very few ranks and large system with many k-points that
        ! could break it. The overlap and Hamiltonian I don't have any buffers for since those will
        ! be constructed and destroyed on-the-fly. ELSI is not smart enough to keep the overlap in
        ! single-proc mode anyway, and the time for Fourier transform is negligible.
        do i=1,KS%n_kpoint
            allocate(KS%kpoint(i)%eigenvalue(KS%n_basis,KS%n_spin))
            allocate(KS%kpoint(i)%occupation(KS%n_basis,KS%n_spin))
            KS%kpoint(i)%eigenvalue=0.0_r8
            KS%kpoint(i)%occupation=0.0_r8
            select type(k=>KS%kpoint(i))
            type is(rl_kspace_eigenproblem_kpoint_real)
                allocate(k%eigenvector( KS%n_basis,KS%n_basis,KS%n_spin ))
                allocate(k%densitymatrix( KS%n_basis,KS%n_basis,KS%n_spin ))
                k%eigenvector=0.0_r8
                k%densitymatrix=0.0_r8
            type is(rl_kspace_eigenproblem_kpoint_complex)
                allocate(k%eigenvector( KS%n_basis,KS%n_basis,KS%n_spin ))
                allocate(k%densitymatrix( KS%n_basis,KS%n_basis,KS%n_spin ))
                k%eigenvector=0.0_r8
                k%densitymatrix=0.0_r8
            end select
        enddo

        ! And now I think we are done, at least for single-proc things!
    !end block sng
    type is(rl_kspace_eigenproblem_multiproc)
    ! Here I set things up for ELSI multi-proc mode, which is pretty much the default
    ! in most production runs.
    !mlt: block

        nprob=kmesh%n_irr_kpoint*KS%n_spin
        call mem%allocate(probind,[2,nprob],persistent=.false.,scalable=.false.)
        call mem%allocate(di,nprob,persistent=.false.,scalable=.false.)
        probind=0
        di=0

        ! So, each problem is defined by a k-point index and a spin index.
        l=0
        do i=1,n_spin
        do j=1,kmesh%n_irr_kpoint
            l=l+1
            probind(:,l)=[j,i]
        enddo
        enddo

        ! Count ranks per problem
        do i=1,nprob
        do j=1,mw%n
            if ( mod(j,nprob) .eq. i-1) di(i)=di(i)+1
        enddo
        enddo

        ! Start assigning ranks to groups. I do it sequentially for some reason,
        ! i.e. that ranks 0-3 gets the first problem, 4-7 the next and so on.
        l=0
        do i=1,nprob
        do j=1,di(i)
            if ( mw%r .eq. l ) then
                KS%group_index=i
                KS%index_irreducible_kpoint=probind(1,i)
                KS%index_spin=probind(2,i)
            endif
            l=l+1
        enddo
        enddo
        ! Keep track of how many groups there are, in total.
        KS%n_group_global=nprob

        ! Now we are ready to split the communicator into smaller chunks!
        ! One chunk per group of ranks that deal with a single eigenproblem.
        call mw%split(KS%ml,KS%group_index)

        ! It is not enough with just the MPI communicater, we need the guy that
        ! handles BLACS as well. If I got it right, the BLACS context should be
        ! aligned with the group communicator. Think this is the case.
        call KS%bl%init(KS%ml)

        ! Make space for distributed matrix storage! First decide
        ! wether it is real or complex.
        if ( complex_eigenvectors ) then
            allocate(rl_kspace_eigenproblem_prob_complex::KS%buf)
        else
            allocate(rl_kspace_eigenproblem_prob_real::KS%buf)
        endif
        ! We have to initialize the block-cyclic storage thingy. For that I need a blocksize.
        KS%buf%blocksize=KS%bl%get_blocksize(KS%n_basis,KS%n_basis)

        ! Now we can initialize ELSI!
        call elsi_init(KS%eh,&
            solver=solver_ELPA,&
            parallel_mode=1,&
            matrix_format=0,&
            n_basis=KS%n_basis,&
            n_electron=KS%n_electron,&
            n_state=KS%n_state)

        ! Sets the local MPI communicator
        call elsi_set_mpi(KS%eh,KS%ml%comm)
        ! And the global MPI communicator
        call elsi_set_mpi_global(KS%eh,mw%comm)

        ! tell ELSI how I distribute matrices
        call elsi_set_blacs(KS%eh,KS%bl%icontxt,KS%buf%blocksize)

        ! gotta tell ELSI how many k-points, spin channels, weight
        ! and stuff like that. Seems elsi gets confused otherwise.
        call elsi_set_kpoint(KS%eh,kmesh%n_irr_kpoint,KS%index_irreducible_kpoint,kmesh%ip(KS%index_irreducible_kpoint)%integration_weight)
        call elsi_set_spin(KS%eh,KS%n_spin,KS%index_spin)

        ! Make space for solution buffers
        KS%buf%n_row_local=KS%bl%count_row_local(KS%n_basis,KS%buf%blocksize)
        KS%buf%n_col_local=KS%bl%count_col_local(KS%n_basis,KS%buf%blocksize)
        ! Small sanity test. Should be imposssible to trigger.
        if ( KS%buf%n_row_local*KS%buf%n_col_local .eq. 0 ) then
            call rl_stop_gracefully(['Empty buffer on one or more ranks, that is not good'],rl_exitcode_mpi,mw%comm)
        endif
        allocate(KS%buf%eigenvalue(KS%n_basis))
        KS%buf%eigenvalue=0.0_r8
        select type(b=>KS%buf)
        type is(rl_kspace_eigenproblem_prob_real)
            allocate(b%hamiltonian  (b%n_row_local,b%n_col_local))
            allocate(b%overlap      (b%n_row_local,b%n_col_local))
            allocate(b%eigenvector  (b%n_row_local,b%n_col_local))
            allocate(b%densitymatrix(b%n_row_local,b%n_col_local))
            b%hamiltonian=  0.0_r8
            b%overlap=      0.0_r8
            b%eigenvector=  0.0_r8
            b%densitymatrix=0.0_r8
        type is(rl_kspace_eigenproblem_prob_complex)
            allocate(b%hamiltonian  (b%n_row_local,b%n_col_local))
            allocate(b%overlap      (b%n_row_local,b%n_col_local))
            allocate(b%eigenvector  (b%n_row_local,b%n_col_local))
            allocate(b%densitymatrix(b%n_row_local,b%n_col_local))
            b%hamiltonian=  0.0_r8
            b%overlap=      0.0_r8
            b%eigenvector=  0.0_r8
            b%densitymatrix=0.0_r8
        end select

        call mem%deallocate(probind,persistent=.false.,scalable=.false.)
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    !end block mlt
    end select

    ! Set the things that don't depend on parallelism
    ! Also likely a good place to look for conflicting settings.
    !setelsi: block
        call elsi_set_sing_check(KS%eh,1) ! Check overlap singularity. Not sure what 1 means.
        !call elsi_set_sing_stop(KS%eh,1)  ! Always stop if overlap is singular, unless override_illconditioning
        call elsi_set_sing_tol(KS%eh,basis_threshold) ! Singularity threshold
        call elsi_set_zero_def(KS%eh,zero_threshold)  ! Numerical zero threshold
        if ( mw%talk ) then
            call elsi_set_output_log(KS%eh,0)
            !call elsi_set_output(KS%eh,elsi_out_level)
            call elsi_set_output(KS%eh,0)
        else
            call elsi_set_output_log(KS%eh,0)
            call elsi_set_output(KS%eh,0)
        endif
        call elsi_set_write_unit(KS%eh,rl_iou) ! not sure what unit to use here. zero for now.
        ! Versioning information. Will not mess with.
        !call write_bare_aims_uuid(uuid)
        !call elsi_set_uuid(eh_scf,trim(uuid))
        ! Broadening scheme and width
        call elsi_set_mu_broaden_scheme(KS%eh,occupation_type)
        if(occupation_type == 4) then ! Cubic polynomial
            call elsi_set_mu_broaden_scheme(KS%eh,3)
        elseif(occupation_type == 5) then ! Cold
            call elsi_set_mu_broaden_scheme(KS%eh,4)
        endif
        call elsi_set_mu_broaden_width(KS%eh,occupation_width)
        call elsi_set_mu_tol(KS%eh,occupation_acc)
        !if(occupation_type == 2) then ! Methfessel-Paxton
        !    call elsi_set_mu_mp_order(eh_scf,n_methfessel_paxton)
        !endif

        ! Solver settings
        select case(elsi_solver)
        case(1)
            if(solver_method /= 0) then
                call elsi_set_elpa_solver(KS%eh,solver_method)
            endif
            ! if(use_gpu_elpa) then
            !   call elsi_set_elpa_gpu(eh_scf,1)
            ! endif
            ! if(use_gpu_kernels_elpa) then
            !   call elsi_set_elpa_gpu_kernels(eh_scf,1)
            ! endif
            call elsi_set_elpa_n_single(KS%eh,elsi_elpa_n_single)
        case default
            call rl_stop_gracefully(['Have not accomodated all ELSI things yet'],rl_exitcode_param,mw%comm)
        end select

    !end block setelsi

    ! Report some things to where we should write things
    !report: block
        if ( mw%talk ) then
            write(*,*) '... finished setting up space for solution of KS eigenproblem (',tochar(mpi_wtime()-timer),'s)'
        endif
    !end block report

    ! Check that I did not do anything stupid
    if ( mem%temporary_scalable .ne. 0 ) call rl_stop_gracefully(['temporary scalable memory not cleared'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 ) call rl_stop_gracefully(['temporary nonscalable memory not cleared'],rl_exitcode_memory,mw%comm)

end subroutine

end module
