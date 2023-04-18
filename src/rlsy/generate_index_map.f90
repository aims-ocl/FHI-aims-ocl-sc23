module generate_index_map
!!
!! So I changed the purpose of this module to provide interfaces that are easy to make
!! stubs out of. It's kinda dumb to just send things along to another routine, but
!! that's what you get when you are not allowed to use precompilers.
!!
use rlsy_constants, only: r8,rl_iou
use rlsy_interface, only: rlsy_handle
use rlsy_mpi_helper, only: rl_stop_gracefully
use rlsy_calculate_density, only: rl_get_density_from_densitymatrix
use rlsy_calculate_overlap, only: rl_calculate_overlapmatrix
use rlsy_solve_KS_equations, only: rl_solve_KS
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc,&
    rl_kspace_eigenproblem_prob,rl_kspace_eigenproblem_prob_real,&
    rl_kspace_eigenproblem_prob_complex,rl_kspace_eigenproblem_kpoint_real,&
    rl_kspace_eigenproblem_kpoint_complex
use rlsy_hartree_potential, only: rl_calculate_multipole_expansion
use grids, only: batch_of_points

implicit none
private

public :: rl_generate_index_map
public :: rl_collect_overlap
public :: rl_collect_hamiltonian
public :: rl_collect_densitymatrix
public :: rl_solve_KS_equations
public :: rl_multipole_expansion
public :: rl_calculate_density
public :: rl_calculate_overlap
public :: rl_check_if_reinit_needed
public :: rl_inject_density
public :: rl_inject_KS_solution
public :: rl_inject_densitymatrix

contains

! Note: There is a good reason why I have separate routines for
! collecting the overlap, hamiltonian, whatever. It's all copy-pasted
! boilerplate code, but there is no place in the AIMS scf cycle that I
! can be sure that all of these quantities are defined at the same time
! so I have to grab things at different points.

!> check if a re-initialization of the symmetry handle is needed. Will have to add more things as we progress.
subroutine rl_check_if_reinit_needed(rlsy_h,frac_coords,lattice_vector,need_reinit)
    type(rlsy_handle), intent(inout) :: rlsy_h
    real(r8), dimension(:,:), intent(in) :: frac_coords
    real(r8), dimension(3,3), intent(in) :: lattice_vector
    logical, intent(out) :: need_reinit

    call rlsy_h%check_handle_rebuild(frac_coords,lattice_vector,need_reinit)
end subroutine

!> This routine collects the Overlap matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_overlap(rlsy_h,overlap_matrix)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real(r8), dimension(:), intent(in) :: overlap_matrix

    !call rlsy_h%rmtx%collect_overlap(rlsy_h%basis,overlap_matrix)
    call rlsy_h%rmtx%collect_overlap(rlsy_h%basis,rlsy_h%structure,rlsy_h%KS,rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,overlap_matrix)
end subroutine

!> This routine collects the Overlap matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_densitymatrix(rlsy_h,densitymatrix)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real(r8), dimension(:), intent(in) :: densitymatrix

    call rlsy_h%rmtx%collect_densitymatrix(rlsy_h%basis,densitymatrix)
end subroutine

!> This routine collects the Hamiltonian matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_hamiltonian(rlsy_h,hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real(r8), dimension(:,:), intent(in) :: hamiltonian

    !call rlsy_h%rmtx%collect_hamiltonian(rlsy_h%basis,hamiltonian)
    call rlsy_h%rmtx%collect_hamiltonian(rlsy_h%basis,rlsy_h%structure,rlsy_h%KS,rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,hamiltonian)
end subroutine

!> This routine injects the Hamiltonian back into AIMS. It assumes that the index is already generated.
subroutine rl_inject_hamiltonian(rlsy_h,hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real(r8), dimension(:,:), intent(inout) :: hamiltonian

    call rlsy_h%rmtx%inject_hamiltonian(rlsy_h%basis,rlsy_h%structure,rlsy_h%KS,rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,hamiltonian)
end subroutine

!> This routine injects the Hamiltonian back into AIMS. It assumes that the index is already generated.
subroutine rl_inject_densitymatrix(rlsy_h,density_matrix_sparse,ispin)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> density matrix
    real(r8), dimension(:), intent(inout) :: density_matrix_sparse
    !> spin index
    integer, intent(in) :: ispin

    !subroutine inject_densitymatrix(rmtx,basis,p,KS,sym,mw,mem,density_matrix_sparse,ispin)
    call rlsy_h%rmtx%inject_densitymatrix(rlsy_h%basis,rlsy_h%structure,rlsy_h%KS,rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,density_matrix_sparse,ispin)
end subroutine

!> This is very annoying, but I have to re-index the Hamiltonian into something useful.
subroutine rl_generate_index_map(rlsy_h,&
    coords_center, center_to_cell, center_to_atom, &
    index_hamiltonian, column_index_hamiltonian, cbasis_to_basis, &
    cbasis_to_center, centers_basis_integrals, frac_coords, lattice_vector,&
    position_in_hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> things from AIMS such that I can figure out the indexing:
    !> AIMS coords_centers, Where is center i
    real(r8), intent(in), dimension(:,:) :: coords_center
    !> AIMS center_to_cell. What vector moves center i back to the unit cell
    integer, intent(in), dimension(:) :: center_to_cell
    !> AIMS center_to_atom, what index in the unit cell is center i?
    integer, intent(in), dimension(:) :: center_to_atom
    !> AIMS index Hamiltonian
    integer, intent(in), dimension(:,:,:) :: index_hamiltonian
    !> AIMS column index Hamiltonian
    integer, intent(in), dimension(:) :: column_index_hamiltonian
    !> AIMS indices that map the all the basis functions to something I can understand.
    integer, intent(in), dimension(:) :: cbasis_to_basis
    integer, intent(in), dimension(:) :: cbasis_to_center
    integer, intent(in), dimension(:) :: centers_basis_integrals
    !> AIMS coordinates of atoms
    real(r8), intent(in), dimension(:,:) :: frac_coords
    real(r8), intent(in), dimension(3,3) :: lattice_vector
    !> AIMS position_in_hamiltonian
    integer, dimension(:,:), intent(in) :: position_in_hamiltonian

    ! It is a little stupid to just have a wrapper subroutine that calls another
    ! one. But I don't see another way since precompilers are apparently not a thing,
    ! and emulating a precompiler with cmake is apparently a much better solution.
    call rlsy_h%rmtx%build_index(rlsy_h%structure,rlsy_h%basis,rlsy_h%mw,rlsy_h%mem,&
            coords_center, center_to_cell, center_to_atom, &
            index_hamiltonian, column_index_hamiltonian, cbasis_to_basis, &
            cbasis_to_center, centers_basis_integrals, frac_coords, lattice_vector,&
            position_in_hamiltonian,&
            rlsy_h%verbosity)
end subroutine

!> This routine injects the Hamiltonian back into AIMS. It assumes that the index is already generated.
subroutine rl_solve_KS_equations(rlsy_h)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Send things along, same procedure as always!
    call rl_solve_KS(rlsy_h%rmtx,rlsy_h%KS,rlsy_h%structure,rlsy_h%basis,rlsy_h%kmesh,&
        rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,rlsy_h%verbosity+1,rlsy_h%timer%eigenvalue)


end subroutine

!> evaluates the density
subroutine rl_calculate_density(rlsy_h)
    !> MPI helper
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Just send things along:
    call rl_get_density_from_densitymatrix(rlsy_h%density,rlsy_h%grid,rlsy_h%rmtx,rlsy_h%KS,rlsy_h%structure,&
        rlsy_h%basis,rlsy_h%ec,rlsy_h%spacegroup,rlsy_h%mw,rlsy_h%mem,rlsy_h%verbosity,rlsy_h%timer%density)
end subroutine

!> evaluates the multipole expansion of the density
subroutine rl_multipole_expansion(rlsy_h)
    !> MPI helper
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Just send things along:
    call rl_calculate_multipole_expansion(rlsy_h%density,rlsy_h%grid,rlsy_h%structure,rlsy_h%spacegroup,&
        rlsy_h%basis,rlsy_h%free,rlsy_h%ec,rlsy_h%mw,rlsy_h%mem,rlsy_h%verbosity,rlsy_h%timer%multipole)
end subroutine


!> evaluates the overlap
subroutine rl_calculate_overlap(rlsy_h)
    !> MPI helper
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Just send things along:
    call rl_calculate_overlapmatrix(rlsy_h%grid,rlsy_h%rmtx,rlsy_h%KS,rlsy_h%structure,&
        rlsy_h%basis,rlsy_h%ec,rlsy_h%mw,rlsy_h%mem,rlsy_h%verbosity,rlsy_h%timer%overlap)
end subroutine

!> replace the AIMS density and gradient of density with the one I have calculated.
subroutine rl_inject_density(rlsy_h,&
    rho,rho_gradient,partition_tab,&
    n_aims_batches,aims_batches)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Densities in different ways from AIMS
    real(r8), dimension(:,:), intent(inout) :: rho
    real(r8), dimension(:,:,:), intent(inout) :: rho_gradient
    !> grid things from AIMS
    real(r8), dimension(:), intent(in) :: partition_tab
    !> Normal batches from AIMS
    integer, intent(in) :: n_aims_batches
    type(batch_of_points), dimension(:), intent(in) :: aims_batches

    ! This is a magical parameter I should choose in a smarter
    ! way. Irrelevant now, but could become important later.
    integer, parameter :: n_pts_per_iter=100000
    logical :: collect_gradient
    ! buffers for full values
    real(r8), dimension(:,:,:), allocatable :: buf_rho_grad
    real(r8), dimension(:,:), allocatable :: buf_rho

    real(r8), dimension(3) :: v0
    real(r8) :: f0,f1,f2
    integer :: niter,iter,ifull
    integer :: i,ib,ip,jp,iglob,iglob_offset
    integer :: iatom,irad,iang,ispc,ispin,iop

    ! Check wether I have a density gradient to collect
    ! I wonder if this is a safe test. Possibly.
    if ( size(rho_gradient) .eq. 3*size(rho) ) then
        collect_gradient=.true.
    else
        collect_gradient=.false.
    endif

    ! Figure out how many times we have to iterate to fill out all the values.
    niter=0
    do
        if ( niter*n_pts_per_iter .ge. rlsy_h%grid%n_point_global ) then
            exit
        else
            niter=niter+1
        endif
    enddo

    ! Space for buffers
    allocate(buf_rho_grad(3,n_pts_per_iter,rlsy_h%KS%n_spin))
    allocate(buf_rho(n_pts_per_iter,rlsy_h%KS%n_spin))
    buf_rho_grad=0.0_r8
    buf_rho=0.0_r8

    ! I iterate over the whole rlsy_h%grid many times, to avoid having to allocate arrays
    ! That are as large as the full number of points on the rlsy_h%grid, since that would
    ! end badly very quickly. So only a subset of points are accumulated each
    ! iteration, and summed up at the end.
    iterloop1: do iter=1,niter
        ! Which global indices do I care about this iteration?
        iglob_offset=(iter-1)*n_pts_per_iter
        ! Reset buffers and counters
        buf_rho=0.0_r8
        buf_rho_grad=0.0_r8

        ! Go through my batches and fill out the relevant entries for this iteration.
        do ib=1,rlsy_h%grid%n_irr_batch
        do ip=1,rlsy_h%grid%irr_batch(ib)%n_point
            jp=rlsy_h%grid%irr_batch(ib)%semilocal_irr_offset+ip
            do i=1,rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
                iatom=rlsy_h%grid%irr_batch(ib)%unfold_atom(i,ip)
                irad=rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i,ip)
                iang=rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i,ip)
                iop=rlsy_h%grid%irr_batch(ib)%unfold_operation(i,ip)
                iglob=rlsy_h%grid%idx%global_index(iatom,irad,iang)-iglob_offset
                ! Check if this point is relevant
                if ( iglob .gt. n_pts_per_iter .or. iglob .le. 0 ) cycle
                ! Now collect scalar things in the right place.
                do ispin=1,rlsy_h%KS%n_spin
                    buf_rho(iglob,ispin)=rlsy_h%density%irr_rho(jp,ispin)
                enddo
                ! Vectorial things, a little worse
                do ispin=1,rlsy_h%KS%n_spin
                    v0=rlsy_h%density%irr_grad_rho(:,jp,ispin)
                    buf_rho_grad(:,iglob,ispin)=matmul(rlsy_h%spacegroup%op(iop)%m,v0)
                enddo
            enddo
        enddo
        enddo
        ! Sync across ranks
        call rlsy_h%mw%allreduce('sum',buf_rho)
        call rlsy_h%mw%allreduce('sum',buf_rho_grad)
        ! Now set these numbers in the AIMS arrays
        ifull=0
        do ib=1,n_aims_batches
        do ip=1,aims_batches(ib)%size
            ! keep track of full index in the aims arrays
            ifull=ifull+1
            ! get global index of this integration point
            iatom=aims_batches(ib)%points(ip)%index_atom
            ispc=rlsy_h%structure%species(iatom)
            irad=aims_batches(ib)%points(ip)%index_radial
            iang=aims_batches(ib)%points(ip)%index_angular
            iglob=rlsy_h%grid%idx%global_index(iatom,irad,iang)-iglob_offset
            ! skip if outside range, get it on the next iteration instead.
            if ( iglob .gt. n_pts_per_iter .or. iglob .le. 0 ) cycle
            ! Store some things
            do ispin=1,rlsy_h%KS%n_spin
                rho(ispin,ifull)=buf_rho(iglob,ispin)
                if ( collect_gradient ) then
                    rho_gradient(:,ispin,ifull)=buf_rho_grad(:,iglob,ispin)
                endif
            enddo
        enddo
        enddo
    enddo iterloop1

    ! Might be sensible to check wether it is properly normalized.
    f0=0.0_r8
    ifull=0
    do ib=1,n_aims_batches
    do ip=1,aims_batches(ib)%size
        ifull=ifull+1
        f0=f0+sum(rho(:,ifull))*partition_tab(ifull)
    enddo
    enddo
    call rlsy_h%mw%allreduce('sum',f0)

    f1=1.0_r8-rlsy_h%KS%n_electron/f0
    if ( abs(f1) .gt. 1E-14_r8 ) then
        ! probably worth it to normalize
        f2=rlsy_h%KS%n_electron/f0
        rho=rho*f2
        rho_gradient=rho_gradient*f0
    endif
end subroutine

!> replace the current solution to the KS equations with the one I have calculated
subroutine rl_inject_KS_solution(rlsy_h,KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,chemical_potential,occ_numbers,k_point_list)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Things from AIMS:
    real(r8), dimension(:,:,:), intent(inout) :: KS_eigenvalue
    real(r8), dimension(:,:,:,:), intent(inout) :: KS_eigenvector
    complex(r8), dimension(:,:,:,:), intent(inout) :: KS_eigenvector_complex
    real(r8), intent(inout) :: chemical_potential
    real(r8), dimension(:,:,:), intent(inout) :: occ_numbers
    real(r8), dimension(:,:), intent(in) :: k_point_list

    integer :: ikp,jkp,kkp
    integer :: i,j,k,l


    ! This is really quite annoying at the moment. Everything can be
    ! distributed in a number of quite arbitrary ways. Much frustration,
    ! and a lot of indexing that needs to be sorted out.

    if ( rlsy_h%mw%talk ) then
        write(*,*) 'INJECTING KS SOLUTION'
        write(*,*) '             shp eigenvalue:',shape(KS_eigenvalue)
        write(*,*) '       shp eigenvector_real:',shape(KS_eigenvector)
        write(*,*) '    shp eigenvector_complex:',shape(KS_eigenvector_complex)
    endif

    ! First the KS eigenvalues. Man this is annoying. Hmm.
    KS_eigenvalue=0.0_r8
    select type(ks=>rlsy_h%KS)
    type is(rl_kspace_eigenproblem_singleproc)
        do i=1,ks%n_kpoint
        enddo
    type is(rl_kspace_eigenproblem_multiproc)
        write(*,*) 'CONTINUE HERE INJECT KS'
        call rlsy_h%mw%destroy()
        stop
    end select

    ! First thing we need to do is match the AIMS k-points with mine.


    write(*,*) 'CONTINUE HERE INJECT KS'
    call rlsy_h%mw%destroy()
    stop
end subroutine

end module
