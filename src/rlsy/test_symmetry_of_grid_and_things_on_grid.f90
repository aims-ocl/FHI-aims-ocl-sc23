module test_symmetry_of_grid_and_things_on_grid
!!
!! Collection of routines to measury how symmetric, or not, properties in AIMS are.
!! This is intended for debugging/diagnostics. For production runs doing these tests
!! serves no purpose whatsoever, so don't do that.
!!
use rlsy_constants, only: r8,rl_iou,rl_exitcode_symmetry
use rlsy_helpers, only: rl_clean_fractional_coordinates,norm2
use rlsy_interface, only: rlsy_handle
use rlsy_mpi_helper, only: rl_stop_gracefully

! Normal AIMS modules
use grids, only: batch_of_points
use constants, only: hartree

implicit none
private

public :: rl_test_symmetry_of_grid_and_things_on_grid
public :: rl_test_symmetry_of_density
public :: rl_track_progress

contains

!> progress tracker, for debugging/understanding. More things will be added as needed.
subroutine rl_track_progress(rlsy_h,&
    number_of_loops,total_energy,previous_total_energy,&
    rho_change,ev_sum,previous_ev_sum,diff_forces,diff_stress,&
    AS_stress_on,forces_on)

    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> things grabbed from inside scf_solver
    integer, intent(in) :: number_of_loops
    real(r8), intent(in) :: total_energy
    real(r8), intent(in) :: previous_total_energy
    real(r8), dimension(:), intent(in) :: rho_change
    real(r8), intent(in) :: ev_sum
    real(r8), intent(in) :: previous_ev_sum
    real(r8), intent(in) :: diff_forces
    real(r8), intent(in) :: diff_stress
    logical, intent(in) :: AS_stress_on
    logical, intent(in) :: forces_on

    ! Dump to non-standard unit to be caught from something else.
    if ( rlsy_h%run_as_library .and. rlsy_h%verbosity .gt. 0 ) then
        if ( forces_on .and. AS_stress_on ) then
            write(rlsy_h%stdout,"(2X,'it',I4,1X,F24.14,5(2X,E12.5))") &
                   number_of_loops,total_energy,&
                   sum(abs(rho_change)),&
                   (ev_sum-previous_ev_sum)*hartree,&
                   (total_energy - previous_total_energy)*hartree,&
                   diff_forces,&
                   diff_stress
        elseif ( forces_on ) then
        elseif ( AS_stress_on ) then
        else
            write(rlsy_h%stdout,"(2X,'it',I4,1X,F24.14,3(2X,E12.5),2(2X,A12))") &
                   number_of_loops,total_energy,&
                   sum(abs(rho_change)),&
                   (ev_sum-previous_ev_sum)*hartree,&
                   (total_energy - previous_total_energy)*hartree,&
                   'N/A','N/A'
        endif
        flush(rlsy_h%stdout)
    endif
end subroutine

!> tests the symmetry of the grid, and quantities defined on the grid.
!@TODO Add potential
!@TODO Add kinetic energy density and gradient
!@TODO The annoying stub things removes all the clarity from the implementation, that is really annoying.
! Maybe I should double-wrap things like ELSI does, but that is also a little stupid.
subroutine rl_test_symmetry_of_density(rlsy_h,&
    rho,delta_rho,rho_free_superpos,rho_gradient,delta_rho_gradient,&
    weight_tab,partition_tab,&
    r_radial,r_angular,&
    n_aims_batches,aims_batches)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Densities in different ways from AIMS
    real(r8), dimension(:,:), intent(in) :: rho
    real(r8), dimension(:,:), intent(in) :: delta_rho
    real(r8), dimension(:), intent(in) :: rho_free_superpos
    real(r8), dimension(:,:,:), intent(in) :: rho_gradient
    real(r8), dimension(:,:,:), intent(in) :: delta_rho_gradient
    !> grid things from AIMS
    real(r8), dimension(:), intent(in) :: weight_tab
    real(r8), dimension(:), intent(in) :: partition_tab
    !> coordinates of points from AIMS
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    !> Normal batches from AIMS
    integer, intent(in) :: n_aims_batches
    type(batch_of_points), dimension(:), intent(in) :: aims_batches

    ! This is a magical parameter I should choose in a smarter way. Irrelevant now,
    ! But could become important later.
    integer, parameter :: n_pts_per_iter=50000
    real(r8), dimension(3,rlsy_h%structure%n_atom) :: folded_rcart
    real(r8) :: timer,t0,t1
    ! buffers for the resulting statistics
    real(r8), dimension(rlsy_h%density%n_spin,3) :: stat_rho_grad
    real(r8), dimension(rlsy_h%density%n_spin,3) :: stat_rho
    real(r8), dimension(3) :: stat_free_rho
    real(r8), dimension(3) :: stat_part_tab
    logical :: collect_gradient
    !init: block
        integer :: i,l
    !avg: block
        ! buffers communication and storing values
        real(r8), dimension(:,:,:), allocatable :: buf_rlsy_grad,buf_aims_grad
        real(r8), dimension(:,:), allocatable :: buf_rlsy_rho,buf_aims_rho
        real(r8), dimension(:,:), allocatable :: buf_rlsy_coord,buf_aims_coord
        real(r8), dimension(:), allocatable :: buf_rlsy_free,buf_aims_free
        real(r8), dimension(:), allocatable :: buf_rlsy_part,buf_aims_part


        real(r8), dimension(3) :: v0,v1
        real(r8) :: f0,f1,f2
        integer :: niter,iter,ifull,nspin
        !integer :: i,j,ib,ip,jp,iglob,iglob_offset,ispin,ispc
        integer :: j,ib,ip,jp,iglob,iglob_offset,ispin,ispc
        integer :: iatom,irad,iang,iop !,jatom,jrad,jang,jop

    !init: block

        ! Well-adjusted coordinates
        do i=1,rlsy_h%structure%n_atom
            folded_rcart(:,i)=rl_clean_fractional_coordinates(rlsy_h%structure%fractional_coordinate(:,i)+0.5_r8)-0.5_r8
            folded_rcart(:,i)=matmul(rlsy_h%structure%latticevectors,folded_rcart(:,i))
        enddo

        ! Initialize the statistics to nothing
        stat_rho=0.0_r8
        stat_rho_grad=0.0_r8
        stat_free_rho=0.0_r8
        stat_part_tab=0.0_r8

        ! Check wether I have a density gradient to collect
        ! I wonder if this is a safe test.
        if ( size(rho_gradient) .eq. 3*size(rho) ) then
            collect_gradient=.true.
        else
            collect_gradient=.false.
        endif
    !end block init




    ! Average things over the grid
    !avg: block

        ! How many iterations are needed. This should be enough.
        niter=0
        do
            if ( niter*n_pts_per_iter .ge. rlsy_h%grid%n_point_global ) then
                exit
            else
                niter=niter+1
            endif
        end do
        ! How many spin channels
        nspin=rlsy_h%KS%n_spin

        allocate(buf_rlsy_grad(3,n_pts_per_iter,nspin))
        allocate(buf_aims_grad(3,n_pts_per_iter,nspin))
        allocate(buf_rlsy_rho(n_pts_per_iter,nspin))
        allocate(buf_aims_rho(n_pts_per_iter,nspin))
        allocate(buf_rlsy_coord(3,n_pts_per_iter))
        allocate(buf_aims_coord(3,n_pts_per_iter))
        allocate(buf_rlsy_free(n_pts_per_iter))
        allocate(buf_aims_free(n_pts_per_iter))
        allocate(buf_rlsy_part(n_pts_per_iter))
        allocate(buf_aims_part(n_pts_per_iter))


        ! I iterate over the whole rlsy_h%grid many times, to avoid having to allocate arrays
        ! That are as large as the full number of points on the rlsy_h%grid, since that would
        ! end badly very quickly. So only a subset of points are accumulated each
        ! iteration, and summed up at the end.
        iterloop1: do iter=1,niter
            ! Which global indices do I care about this iteration?
            iglob_offset=(iter-1)*n_pts_per_iter
            ! Reset buffers and counters
            buf_rlsy_grad=0.0_r8
            buf_aims_grad=0.0_r8
            buf_rlsy_rho=0.0_r8
            buf_aims_rho=0.0_r8
            buf_rlsy_coord=0.0_r8
            buf_aims_coord=0.0_r8
            buf_rlsy_free=0.0_r8
            buf_aims_free=0.0_r8
            buf_rlsy_part=0.0_r8
            buf_aims_part=0.0_r8

            ! First we collect from the AIMS mesh
            ifull=0
            do ib=1,n_aims_batches
            do ip=1,aims_batches(ib)%size
                ! keep track of full index in the aims arrays
                ifull=ifull+1
                ! get global index of this integration point
                iatom=aims_batches(ib)%points(ip)%index_atom
                !ispc=rlsy_h%structure%species(iatom)
                irad=aims_batches(ib)%points(ip)%index_radial
                iang=aims_batches(ib)%points(ip)%index_angular
                iglob=rlsy_h%grid%idx%global_index(iatom,irad,iang)-iglob_offset
                ! skip if outside range, get it on the next iteration instead.
                if ( iglob .gt. n_pts_per_iter .or. iglob .le. 0 ) cycle
                ! Ok, this point is within the current range. Store some things.
                do ispin=1,nspin
                    buf_aims_rho(iglob,ispin)=buf_aims_rho(iglob,ispin)+rho(ispin,ifull)+delta_rho(ifull,ispin)
                    if ( collect_gradient ) then
                        buf_aims_grad(:,iglob,ispin)=rho_gradient(:,ispin,ifull)+delta_rho_gradient(:,ifull,ispin)
                    else
                        buf_aims_grad(:,iglob,ispin)=0.0_r8
                    endif
                enddo
                buf_aims_free(iglob)=rho_free_superpos(ifull)              ! free atom density
                buf_aims_part(iglob)=partition_tab(ifull)                  ! partition function
                buf_aims_coord(:,iglob)=aims_batches(ib)%points(ip)%coords ! coordinate, later problem.
            enddo
            enddo

            ! Then collect from my mesh
            do ib=1,rlsy_h%grid%n_irr_batch
            do ip=1,rlsy_h%grid%irr_batch(ib)%n_point
                jp=rlsy_h%grid%irr_batch(ib)%semilocal_irr_offset+ip
                do i=1,rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
                    iatom=rlsy_h%grid%irr_batch(ib)%unfold_atom(i,ip)
                    ispc=rlsy_h%structure%species(iatom)
                    irad=rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i,ip)
                    iang=rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i,ip)
                    iop=rlsy_h%grid%irr_batch(ib)%unfold_operation(i,ip)
                    iglob=rlsy_h%grid%idx%global_index(iatom,irad,iang)-iglob_offset
                    if ( iglob .gt. n_pts_per_iter .or. iglob .le. 0 ) cycle
                    ! Now collect scalar things in the right place.
                    do ispin=1,nspin
                        buf_rlsy_rho(iglob,ispin)=rlsy_h%density%irr_rho(jp,ispin)
                        if ( collect_gradient ) then
                            v0=rlsy_h%density%irr_grad_rho(:,jp,ispin)
                            buf_rlsy_grad(:,iglob,ispin)=matmul(rlsy_h%spacegroup%op(iop)%m,v0)
                        else
                            buf_rlsy_grad(:,iglob,ispin)=0.0_r8
                        endif
                    enddo
                    buf_rlsy_free(iglob)=rlsy_h%density%irr_free_rho(jp)
                    f2=rlsy_h%grid%irr_batch(ib)%partition_function(ip)*rlsy_h%grid%irr_batch(ib)%integration_weight(ip)/real(rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip),r8)
                    buf_rlsy_part(iglob)=f2
                    buf_rlsy_coord(:,iglob)=folded_rcart(:,iatom)+r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)
                enddo
            enddo
            enddo

            ! The buffers should be synced
            call rlsy_h%mw%allreduce('sum',buf_rlsy_grad)
            call rlsy_h%mw%allreduce('sum',buf_aims_grad)
            call rlsy_h%mw%allreduce('sum',buf_rlsy_rho)
            call rlsy_h%mw%allreduce('sum',buf_aims_rho)
            call rlsy_h%mw%allreduce('sum',buf_rlsy_coord)
            call rlsy_h%mw%allreduce('sum',buf_aims_coord)
            call rlsy_h%mw%allreduce('sum',buf_rlsy_free)
            call rlsy_h%mw%allreduce('sum',buf_aims_free)
            call rlsy_h%mw%allreduce('sum',buf_rlsy_part)
            call rlsy_h%mw%allreduce('sum',buf_aims_part)

            ! Now I can start comparing:
            if ( rlsy_h%mw%talk ) then
                ispin=1
                do i=1,min(n_pts_per_iter,rlsy_h%grid%n_point_global)

                    ! Don't print all
                    !if ( mod(i,500) .ne. 0 ) cycle
                    if ( buf_aims_part(i) .lt. 1E-45_r8 ) cycle

                    if ( i .eq. 1100 ) then
                    write(*,*) ''
                    write(*,*) 'point',i
                    write(*,*) '  crda:',buf_aims_coord(:,i)
                    write(*,*) '  crdr:',buf_rlsy_coord(:,i)
                    write(*,*) '   rho:',buf_aims_rho(i,ispin),buf_rlsy_rho(i,ispin),buf_aims_rho(i,ispin)-buf_rlsy_rho(i,ispin),buf_aims_rho(i,ispin)/buf_rlsy_rho(i,ispin)
                    write(*,*) '  free:',buf_aims_free(i),buf_rlsy_free(i),buf_aims_free(i)-buf_rlsy_free(i)
                    write(*,*) '  part:',buf_aims_part(i),buf_rlsy_part(i),buf_aims_part(i)-buf_rlsy_part(i)
                    write(*,*) ' grad a:',buf_aims_grad(:,i,ispin)
                    write(*,*) ' grad r:',buf_rlsy_grad(:,i,ispin)
                    endif
                enddo
            endif
        enddo iterloop1
    !end block avg
end subroutine


!> tests the symmetry of the grid, and quantities defined on the grid.
!@TODO Add potential
!@TODO Add kinetic energy density and gradient
!@TODO The annoying stub things removes all the clarity from the implementation, that is really annoying.
! Maybe I should double-wrap things like ELSI does, but that is also a little stupid.
subroutine rl_test_symmetry_of_grid_and_things_on_grid(rlsy_h,&
    rho,delta_rho,rho_free_superpos,rho_gradient,delta_rho_gradient,&
    weight_tab,partition_tab,&
    r_radial,r_angular,&
    n_aims_batches,aims_batches)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Densities in different ways from AIMS
    real(r8), dimension(:,:), intent(in) :: rho
    real(r8), dimension(:,:), intent(in) :: delta_rho
    real(r8), dimension(:), intent(in) :: rho_free_superpos
    real(r8), dimension(:,:,:), intent(in) :: rho_gradient
    real(r8), dimension(:,:,:), intent(in) :: delta_rho_gradient
    !> grid things from AIMS
    real(r8), dimension(:), intent(in) :: weight_tab
    real(r8), dimension(:), intent(in) :: partition_tab
    !> point things from AIMS
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    !> Normal batches from AIMS
    integer, intent(in) :: n_aims_batches
    type(batch_of_points), dimension(:), intent(in) :: aims_batches

    ! This is a magical parameter I should choose in a smarter way. Irrelevant now,
    ! But could become important later.
    integer, parameter :: n_pts_per_iter=50000
    real(r8), dimension(3,rlsy_h%structure%n_atom) :: folded_rcart
    real(r8) :: timer,t0,t1
    ! buffers for the resulting statistics
    real(r8), dimension(rlsy_h%density%n_spin,3) :: stat_rho_grad
    real(r8), dimension(rlsy_h%density%n_spin,3) :: stat_rho
    real(r8), dimension(3) :: stat_free_rho
    real(r8), dimension(3) :: stat_part_tab
    logical :: collect_gradient
    !init: block
        integer :: i,l
    !avg: block
        real(r8), dimension(:,:,:), allocatable :: buf_rho_grad
        real(r8), dimension(:,:), allocatable :: buf_rho
        real(r8), dimension(:), allocatable :: buf_free_rho,buf_part_tab
        ! buffers for irreducible values
        real(r8), dimension(:,:,:), allocatable :: ibuf_rho_grad
        real(r8), dimension(:,:), allocatable :: ibuf_rho
        real(r8), dimension(:), allocatable :: ibuf_free_rho,ibuf_part_tab
        integer, dimension(:), allocatable :: ictr

        real(r8), dimension(3) :: v0,v1
        real(r8) :: f0,f1,f2
        integer :: niter,iter,ifull
        !integer :: i,j,ib,ip,jp,iglob,jglob,iglob_offset
        integer :: j,ib,ip,jp,iglob,jglob,iglob_offset
        integer :: iatom,irad,iang,ispc,jatom,jrad,jang,jop
    !report: block
        real(r8) :: integrated_volume,integrated_charge

    !init: block

        ! Well-adjusted coordinates
        do i=1,rlsy_h%structure%n_atom
            folded_rcart(:,i)=rl_clean_fractional_coordinates(rlsy_h%structure%fractional_coordinate(:,i)+0.5_r8)-0.5_r8
            folded_rcart(:,i)=matmul(rlsy_h%structure%latticevectors,folded_rcart(:,i))
        enddo

        ! Initialize the statistics to nothing
        stat_rho=0.0_r8
        stat_rho_grad=0.0_r8
        stat_free_rho=0.0_r8
        stat_part_tab=0.0_r8

        ! Check wether I have a density gradient to collect
        ! I wonder if this is a safe test.
        if ( size(rho_gradient) .eq. 3*size(rho) ) then
            collect_gradient=.true.
        else
            collect_gradient=.false.
        endif
    !end block init

    ! Average things over the grid
    !avg: block
        ! buffers for full values

        ! How many iterations are needed. This should be enough.
        niter=ceiling( real(rlsy_h%grid%n_point_global,r8)/real(n_pts_per_iter,r8) )

        ! Space for buffers
        allocate(buf_rho_grad(3,n_pts_per_iter,rlsy_h%density%n_spin))
        allocate(buf_rho(n_pts_per_iter,rlsy_h%density%n_spin))
        allocate(buf_free_rho(n_pts_per_iter))
        allocate(buf_part_tab(n_pts_per_iter))
        buf_rho_grad=0.0_r8
        buf_rho=0.0_r8
        buf_free_rho=0.0_r8
        buf_part_tab=0.0_r8

        allocate(ibuf_rho_grad(3,rlsy_h%density%n_semilocal_point,rlsy_h%density%n_spin))
        allocate(ibuf_rho(rlsy_h%density%n_semilocal_point,rlsy_h%density%n_spin))
        allocate(ibuf_free_rho(rlsy_h%density%n_semilocal_point))
        allocate(ibuf_part_tab(rlsy_h%density%n_semilocal_point))
        allocate(ictr(rlsy_h%density%n_semilocal_point))
        ibuf_rho_grad=0.0_r8
        ibuf_rho=0.0_r8
        ibuf_free_rho=0.0_r8
        ibuf_part_tab=0.0_r8
        ictr=0

        ! I iterate over the whole rlsy_h%grid many times, to avoid having to allocate arrays
        ! That are as large as the full number of points on the rlsy_h%grid, since that would
        ! end badly very quickly. So only a subset of points are accumulated each
        ! iteration, and summed up at the end.
        iterloop1: do iter=1,niter
            ! Which global indices do I care about this iteration?
            iglob_offset=(iter-1)*n_pts_per_iter
            ! Reset buffers and counters
            buf_rho=0.0_r8
            buf_free_rho=0.0_r8
            buf_part_tab=0.0_r8
            buf_rho_grad=0.0_r8
            ifull=0
            ! Start collecting from the full rlsy_h%grid
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
                ! Ok, this point is within the current range. Store some things.
                do i=1,rlsy_h%density%n_spin
                    buf_rho(iglob,i)=rho(i,ifull)+delta_rho(ifull,i)
                    if ( collect_gradient ) then
                        buf_rho_grad(:,iglob,i)=rho_gradient(:,i,ifull)+delta_rho_gradient(:,ifull,i)
                    else
                        buf_rho_grad(:,iglob,i)=0.0_r8
                    endif
                enddo
                buf_free_rho(iglob)=rho_free_superpos(ifull)   ! free atom rlsy_h%density
                buf_part_tab(iglob)=partition_tab(ifull)       ! partition function
            enddo
            enddo
            ! The buffers should be synced
            call rlsy_h%mw%allreduce('sum',buf_rho)
            call rlsy_h%mw%allreduce('sum',buf_rho_grad)
            call rlsy_h%mw%allreduce('sum',buf_free_rho)
            call rlsy_h%mw%allreduce('sum',buf_part_tab)
            ! Now go through the irreducible batches and start building averages
            do ib=1,rlsy_h%grid%n_irr_batch
            do ip=1,rlsy_h%grid%irr_batch(ib)%n_point
                jp=rlsy_h%grid%irr_batch(ib)%semilocal_irr_offset+ip
                do i=1,rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
                    jatom=rlsy_h%grid%irr_batch(ib)%unfold_atom(i,ip)
                    jrad=rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i,ip)
                    jang=rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i,ip)
                    jop=rlsy_h%grid%irr_batch(ib)%unfold_operation(i,ip)
                    jglob=rlsy_h%grid%idx%global_index(jatom,jrad,jang)-iglob_offset
                    if ( jglob .gt. n_pts_per_iter .or. jglob .le. 0 ) cycle
                    ! Now collect scalar things in the right place.
                    do j=1,rlsy_h%density%n_spin
                        ibuf_rho(jp,j)=ibuf_rho(jp,j)+buf_rho(jglob,j)
                    enddo
                    ibuf_free_rho(jp)=ibuf_free_rho(jp)+buf_free_rho(jglob)
                    ibuf_part_tab(jp)=ibuf_part_tab(jp)+buf_part_tab(jglob)
                    ictr(jp)=ictr(jp)+1
                    ! Vectorial things, a little worse
                    do j=1,rlsy_h%density%n_spin
                        v0=buf_rho_grad(:,jglob,j)
                        v0=matmul(rlsy_h%spacegroup%op(jop)%im,v0)
                        ibuf_rho_grad(:,jp,j)=ibuf_rho_grad(:,jp,j)+v0
                    enddo
                enddo
            enddo
            enddo
        enddo iterloop1

        ! Now I have collected things from the full grid to an irreducible representation.
        ! This is a decent place to insert some sanity tests to catch bugs in the symmetry
        ! detection, among other things.
        do ib=1,rlsy_h%grid%n_irr_batch
        do ip=1,rlsy_h%grid%irr_batch(ib)%n_point
            jp=rlsy_h%grid%irr_batch(ib)%semilocal_irr_offset+ip
            ! Test that the multiplicity I have is the one I think I have.
            if ( ictr(jp) .ne. rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip) ) then
                call rl_stop_gracefully(['I did not collect all points'],rl_exitcode_symmetry,rlsy_h%mw%comm)
            endif
            ! If this passes, I can average what I have collected:
            ibuf_rho(jp,:)=ibuf_rho(jp,:)/real(ictr(jp),r8)
            ibuf_free_rho(jp)=ibuf_free_rho(jp)/real(ictr(jp),r8)
            ibuf_part_tab(jp)=ibuf_part_tab(jp)/real(ictr(jp),r8)
            ibuf_rho_grad(:,jp,:)=ibuf_rho_grad(:,jp,:)/real(ictr(jp),r8)
        enddo
        enddo

        ! Now to do the same thing again, but this time I compare the symmetry-averaged
        ! values with the raw values on the grid, and accumulate various integrated things.
        iterloop2: do iter=1,niter
            ! Which global indices do I care about this iteration?
            iglob_offset=(iter-1)*n_pts_per_iter
            ! Reset buffers and counters
            buf_rho=0.0_r8
            buf_free_rho=0.0_r8
            buf_part_tab=0.0_r8
            buf_rho_grad=0.0_r8
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
                do i=1,rlsy_h%density%n_spin
                    buf_rho(iglob,i)=rho(i,ifull)+delta_rho(ifull,i)
                    if ( collect_gradient ) then
                        buf_rho_grad(:,iglob,i)=rho_gradient(:,i,ifull)+delta_rho_gradient(:,ifull,i)
                    else
                        buf_rho_grad(:,iglob,i)=0.0_r8
                    endif
                enddo
                buf_free_rho(iglob)=rho_free_superpos(ifull)   ! free atom rlsy_h%density
                buf_part_tab(iglob)=partition_tab(ifull)       ! partition function
            enddo
            enddo
            ! The buffers should be synced
            call rlsy_h%mw%allreduce('sum',buf_rho)
            call rlsy_h%mw%allreduce('sum',buf_rho_grad)
            call rlsy_h%mw%allreduce('sum',buf_free_rho)
            call rlsy_h%mw%allreduce('sum',buf_part_tab)
            ! Go through and compare with the averaged value
            do ib=1,rlsy_h%grid%n_irr_batch
            do ip=1,rlsy_h%grid%irr_batch(ib)%n_point
                jp=rlsy_h%grid%irr_batch(ib)%semilocal_irr_offset+ip
                do i=1,rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
                    jatom=rlsy_h%grid%irr_batch(ib)%unfold_atom(i,ip)
                    jrad=rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i,ip)
                    jang=rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i,ip)
                    jop=rlsy_h%grid%irr_batch(ib)%unfold_operation(i,ip)
                    jglob=rlsy_h%grid%idx%global_index(jatom,jrad,jang)-iglob_offset
                    if ( jglob .gt. n_pts_per_iter .or. jglob .le. 0 ) cycle
                    ! temp store my partition*weight/symweight
                    f2=rlsy_h%grid%irr_batch(ib)%partition_function(ip)*rlsy_h%grid%irr_batch(ib)%integration_weight(ip)/real(ictr(jp),r8)

                    ! rlsy_h%density
                    do j=1,rlsy_h%density%n_spin
                        f0=abs(buf_rho(jglob,j)-ibuf_rho(jp,j))                ! Absolute deviation from symmetry-average
                        f1=abs(buf_rho(jglob,j)-rlsy_h%density%irr_rho(jp,j) ) ! Deviation from my calculated density
                        stat_rho(j,1)=max(f0,stat_rho(j,1))                ! Largest deviation
                        stat_rho(j,2)=stat_rho(j,2)+f0*buf_part_tab(jglob) ! Integrated deviation
                        stat_rho(j,3)=stat_rho(j,3)+f1*buf_part_tab(jglob) ! Integrated deviation
                    enddo

                    ! free atom rlsy_h%density
                    f0=abs(buf_free_rho(jglob)-ibuf_free_rho(jp))
                    stat_free_rho(1)=max(f0,stat_free_rho(1))
                    stat_free_rho(2)=stat_free_rho(2)+f0*buf_part_tab(jglob)
                    stat_free_rho(3)=stat_free_rho(3)+f0*f2

                    ! partition tab
                    f0=abs(buf_part_tab(jglob)-ibuf_part_tab(jp))
                    stat_part_tab(1)=max(f0,stat_part_tab(1))
                    stat_part_tab(2)=stat_part_tab(2)+f0*buf_part_tab(jglob)
                    stat_part_tab(3)=stat_part_tab(3)+f0*f2

                    ! gradient
                    do j=1,rlsy_h%density%n_spin
                        v0=buf_rho_grad(:,jglob,j)
                        v0=matmul(rlsy_h%spacegroup%op(jop)%im,v0)
                        f0=norm2(ibuf_rho_grad(:,jp,j)-v0)
                        f1=norm2(ibuf_rho_grad(:,jp,j)-rlsy_h%density%irr_grad_rho(:,jp,j))
                        !
                        v1=rlsy_h%density%irr_grad_rho(:,jp,j)
                        !f1=norm2(rlsy_h%density%irr_grad_rho(:,jp,j))
                        stat_rho_grad(j,1)=max(f0,stat_rho_grad(j,1))                           ! Largest deviation
                        stat_rho_grad(j,2)=stat_rho_grad(j,2)+f0*buf_part_tab(jglob)            ! Integrated deviation
                        stat_rho_grad(j,3)=stat_rho_grad(j,3)+f1*buf_part_tab(jglob) ! Integrated deviation my part tab
                    enddo
                enddo
            enddo
            enddo
        enddo iterloop2

        ! Sync up the diagnostics
        call rlsy_h%mw%allreduce('max',stat_rho(:,1))
        call rlsy_h%mw%allreduce('sum',stat_rho(:,2:3))

        call rlsy_h%mw%allreduce('max',stat_free_rho(1))
        call rlsy_h%mw%allreduce('sum',stat_free_rho(2:3))

        call rlsy_h%mw%allreduce('max',stat_part_tab(1))
        call rlsy_h%mw%allreduce('sum',stat_part_tab(2:3))

        call rlsy_h%mw%allreduce('max',stat_rho_grad(:,1))
        call rlsy_h%mw%allreduce('sum',stat_rho_grad(:,2:3))

        ! Some cleanup
        deallocate(buf_rho)
        deallocate(buf_free_rho)
        deallocate(buf_part_tab)
        deallocate(ibuf_rho)
        deallocate(ibuf_free_rho)
        deallocate(ibuf_part_tab)
        deallocate(ictr)
    !end block avg

    ! Report what I found
    !report: block


        integrated_volume=sum(partition_tab)/rlsy_h%structure%volume
        call rlsy_h%mw%allreduce('sum',integrated_volume)
        integrated_charge=sum(partition_tab*(rho(1,:)+delta_rho(:,1)))
        call rlsy_h%mw%allreduce('sum',integrated_charge)

        if ( rlsy_h%mw%talk ) then
            write(rl_iou,*) ''
            write(rl_iou,*) 'Report on rlsy symmetry analysis of the integration grid:'
            write(rl_iou,*) '          Volume:',integrated_volume,' (should be 1 in an ideal world)'
            write(rl_iou,*) '          Charge:',integrated_charge,' (should be N electrons)'
            write(rl_iou,*) 'Electron density'
            write(rl_iou,*) '   max deviation:',stat_rho(:,1),' (should be zero)'
            write(rl_iou,*) '      integrated:',stat_rho(:,2),' (should be zero)'
            write(rl_iou,*) '        wrt mine:',stat_rho(:,3),' (should be zero)'
            write(rl_iou,*) 'Electron density gradient'
            write(rl_iou,*) '   max deviation:',stat_rho_grad(:,1),' (should be zero)'
            write(rl_iou,*) '      integrated:',stat_rho_grad(:,2),' (should be zero)'
            write(rl_iou,*) '        wrt mine:',stat_rho_grad(:,3),' (should be zero)'
            write(rl_iou,*) 'Free atom superposition'
            write(rl_iou,*) '   max deviation:',stat_free_rho(1),' (MUST be zero)'
            write(rl_iou,*) '      integrated:',stat_free_rho(2),' (MUST be zero)'
            write(rl_iou,*) 'Partition function'
            write(rl_iou,*) '   max deviation:',stat_part_tab(1),' (MUST be zero)'
            write(rl_iou,*) '      integrated:',stat_part_tab(2),' (MUST be zero)'
            write(rl_iou,*) ''
        endif
    !end block report
end subroutine

end module
