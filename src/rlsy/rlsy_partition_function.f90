module rlsy_partition_function
!!
!! Handles partitions of unity, poorly. Only for debugging.
!!
use rlsy_constants, only: r8,rl_huge,rl_sqtol,rl_exitcode_param,rl_exitcode_symmetry,rl_pi
use rlsy_helpers, only: rl_clean_fractional_coordinates,rl_chop,norm2
use rlsy_mpi_helper, only: rl_stop_gracefully
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_distancetable, only: rl_distancetable
use rlsy_extended_cluster, only: rl_extended_cluster,rl_extended_cluster_hashedspace

implicit none
private
public :: rl_stratmann_smoother_partition_function

!@TODO probably add more than one kind of partition function.

!> magic stratmann parameter
real(r8), parameter :: stratmann_a=0.64_r8
!> hard-coded onset for Cos smoothing function
real(r8), parameter :: smooth_onset=0.8_r8,smooth_reg=0.2_r8
!> tolerance for what a zero in the weight function is
real(r8), parameter :: wfuntol=1E-10_r8

contains

! Pretty much a copy of the stratmann smoother thingy in AIMS, but
! able to evaluate at some arbitrary point so that I can test symmetry
! and stuff. Also I don't store the distance table, everything is evaluated
! on-the-fly instead, seems to almost be faster and does not use any
! significant amount of memory.
subroutine rl_stratmann_smoother_partition_function(ec,hs,p,iatom,point,cutoff,inner_cutoff_sq,partitionfunction)
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> pre-hashed space that makes distance lookups fast
    type(rl_extended_cluster_hashedspace), intent(inout) :: hs
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> which unit cell atom owns the integration point, i.e. should have it's weight determined.
    integer, intent(in) :: iatom
    !> point to determine partition function for, in absolute Cartesian coordinates
    real(r8), dimension(3), intent(in) :: point
    !> cutoff radius, will only consider contributions this far from iatom.
    real(r8), intent(in) :: cutoff
    !> inner cutoff per atom, squared
    real(r8), dimension(:), intent(in) :: inner_cutoff_sq
    !> resulting partition function for iatom
    real(r8), intent(out) :: partitionfunction
    !> cutoff per atom
    !real(r8), dimension(:), intent(in) :: cutoff_per_atom

    real(r8), parameter :: tinybuf=1E-10_r8
    integer, dimension(hs%max_n_pts) :: didsomething
    real(r8), dimension(hs%max_n_pts) :: weightfunction
    real(r8), dimension(3) :: point_shifted,atom_shifted
    real(r8) :: cutoffsq
    integer :: icluster,n_cluster
    !init: block
        real(r8), dimension(3) :: vl,vr
        real(r8) :: f0,f1
        integer :: i,j
    !partition: block
        !real(r8) :: r_ij,r_ig,r_jg,mu_ij,mu_ji,f0
        real(r8) :: r_ij,r_ig,r_jg,mu_ij,mu_ji
        real(r8) :: fij,fji,gij,gji,sij,sji
        !integer :: i,j

    ! Set up some basic things.
    !init: block

        ! First sanity tests. Check that space is hashed and allocated and so on.
        if ( allocated(hs%box) .eqv. .false. ) then
            call rl_stop_gracefully(['You have to initialize the space-hashing thingy first.'],rl_exitcode_param)
        endif
        if ( hs%cutoff .lt. cutoff ) then
            call rl_stop_gracefully(['Cutoff in hashing not sufficient.'],rl_exitcode_param)
        endif

        ! Squared cutoff for distance checks
        cutoffsq=cutoff**2+1E-10_r8
        ! Initialize partition function to zero
        partitionfunction=0.0_r8

        ! Start transforming the point, and work out some indices. A little
        ! confusing, I admit, but it's all fine. First get the coordinate of iatom,
        ! in the reference frame of the extended cluster. This makes sure we have
        ! iatom such that it is surrounded by a ton of neighbours.
        vr=ec%cartesian_coordinate(:,ec%locate_unit_cell_atom_in_cluster(iatom))
        ! First quick exit, the point might be very close to iatom, in that case
        ! we are done and can return early!
        f0=sqnorm(vr-point)
        if ( f0 .lt. inner_cutoff_sq(iatom) ) then
            ! we know the partition function!
            partitionfunction=1.0_r8
            return
        endif
        if ( f0 .gt. cutoffsq ) then
            ! Point too far away. This should never happen.
            ! Or maybe it should? I dunno. Seems weird, but possibly?
            ! It seems stupid of AIMS to generate radial points that are
            ! really far away, but maybe it happens?
            partitionfunction=0.0_r8
            return
            ! Might put this back in later.
            !call rl_stop_gracefully(['Point outside cutoff, should never happen.'],rl_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! conver the point to fractional coordinates
        vl=matmul(p%inv_latticevectors,point)
        ! shift point with lattice vector to the correct interval, -0.5 <= r < 0.5
        vl=rl_clean_fractional_coordinates(vl+0.5_r8)-0.5_r8

        ! Use the point in fractional correctly shifted coordinates
        ! to fetch the relevant parts of the cluster? The temporary arrays
        ! with distance, index and so on are stored in 'hs'.
        call hs%fetch_relevant_atoms(ec,p,vl,cutoff+tinybuf,n_cluster)
        ! Here is an odd corner case: in case of an atom with no neighbours,
        ! we have to assign the grid point to that atom. That can be detected
        ! via the fact that the cluster only has one atom.
        if ( n_cluster .eq. 1 ) then
            partitionfunction=1.0_r8
            return
        endif

        ! convert point back to Cartesian coordinates
        vl=matmul(p%latticevectors,vl)
        ! subtract original point, this should leave the lattice vector that was used
        ! to shift the point to the right interval
        vl=vl-point
        ! Check that I did that right.
        vl=matmul(p%inv_latticevectors,vl)
        if ( sum(abs(vl-anint(vl))) .gt. rl_sqtol ) then
            call rl_stop_gracefully(['I do not understand vectors'],rl_exitcode_symmetry)
        else
            vl=matmul(p%latticevectors,anint(vl))
        endif
        ! Yup, seems ok. Keep coordinates of shifted point.
        point_shifted=point+vl
        ! Apply the same shift to the atom.
        atom_shifted=vr+vl

        ! Locate the index of the shifted iatom in the cluster:
        icluster=-1
        do i=1,n_cluster
            f1=sqnorm(hs%temp_cartesian_coordinate(:,i)-atom_shifted)
            if ( f1 .lt. rl_sqtol ) then
                ! I have located atom in cluster!
                icluster=i
                exit
            endif
        enddo
        ! Safety check
        if ( icluster .le. 0 ) then
            call rl_stop_gracefully(['I do not fully understand clusters and vectors'],rl_exitcode_symmetry)
        endif

        ! Look for a quick exit, if the point is very close to some other atom we can
        ! exit early, but only if that point is not iatom.
        do i=1,n_cluster
            if ( i .eq. icluster ) cycle
            j=hs%temp_unitcell_index(i)
            if ( hs%temp_dist_sq(i) .lt. inner_cutoff_sq(j) ) then
                ! Quick exit
                partitionfunction=0.0_r8
                return
            endif
        enddo
    !end block init

    ! Begin operation actual operation. If we made it here it means we have to
    ! calculate the actual partition function.
    !partition: block

        weightfunction=1.0_r8
        didsomething=0

        il: do i=1,n_cluster
            r_ig=hs%temp_dist(i)
            ! This should not be a necessary check
            if ( r_ig .gt. cutoff ) cycle il

            jl: do j=i+1,n_cluster
                r_jg=hs%temp_dist(j)
                ! This should also never happen, I think.
                if ( r_jg .gt. cutoff ) cycle jl

                ! Now convert to actual distances.
                r_ij=norm2(hs%temp_cartesian_coordinate(:,i)-hs%temp_cartesian_coordinate(:,j))
                ! And the weird elliptical coordinates.
                mu_ij=(r_ig-r_jg)/r_ij
                mu_ji=(r_jg-r_ig)/r_ij

                ! Stratmann functions
                if ( mu_ij .le. -stratmann_a ) then
                    gij=-1.0_r8
                elseif ( mu_ij .ge. stratmann_a ) then
                    gij=1.0_r8
                else
                    f0=mu_ij/stratmann_a
                    gij=(35*f0 - 35*f0**3 + 21*f0**5 - 5*f0**7)*0.0625_r8
                endif
                if ( mu_ji .le. -stratmann_a ) then
                    gji=-1.0_r8
                elseif ( mu_ji .ge. stratmann_a ) then
                    gji=1.0_r8
                else
                    f0=mu_ji/stratmann_a
                    gji=(35*f0 - 35*f0**3 + 21*f0**5 - 5*f0**7)*0.0625_r8
                endif

                ! Get scaling factors
                if ( r_jg .gt. cutoff ) then
                    fij=0.0_r8
                elseif ( r_jg .lt. smooth_onset*cutoff ) then
                    fij=1.0_r8
                else
                    ! scale to 0-pi
                    fij=rl_pi*(r_jg-smooth_onset*cutoff)/(smooth_reg*cutoff)
                    ! scale to 1-0
                    fij=0.5_r8+cos(fij)*0.5_r8
                endif
                if ( r_ig .gt. cutoff ) then
                    fji=0.0_r8
                elseif ( r_ig .lt. smooth_onset*cutoff ) then
                    fji=1.0_r8
                else
                    ! scale to 0-pi
                    fji=rl_pi*(r_ig-smooth_onset*cutoff)/(smooth_reg*cutoff)
                    ! scale to 1-0
                    fji=0.5_r8+cos(fji)*0.5_r8
                endif

                ! Don't want weird numerics making life annoying.
                gij=min(gij,1.0_r8)
                gij=max(gij,-1.0_r8)
                gji=min(gji,1.0_r8)
                gji=max(gji,-1.0_r8)

                ! Get the cluster function
                sij=0.5_r8*(1-gij)
                sji=0.5_r8*(1-gji)

                ! Scale cluster functions, 'stratmann smoother'
                sij=(1-fij)+fij*sij
                sji=(1-fji)+fji*sji

                sij=max(sij,0.0_r8)
                sij=min(sij,1.0_r8)
                sji=max(sji,0.0_r8)
                sji=min(sji,1.0_r8)

                weightfunction(i)=weightfunction(i)*sij
                weightfunction(j)=weightfunction(j)*sji
                didsomething(i)=1
                didsomething(j)=1
            enddo jl
        enddo il

        weightfunction=weightfunction*didsomething
        !weightfunction=rl_chop(weightfunction,wfuntol)
        weightfunction=weightfunction/sum(weightfunction)
        partitionfunction=weightfunction(icluster)
        !if ( abs(partitionfunction) .lt. wfuntol ) partitionfunction=0.0_r8
        !if ( abs(1-partitionfunction) .lt. wfuntol ) partitionfunction=1.0_r8
    !end block partition
end subroutine

! We do a lot of squared distance checks. I put the function here to ensure it will get inlined.
pure function sqnorm(v) result(nrm)
    real(r8), dimension(3), intent(in) :: v
    real(r8) :: nrm
    nrm=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
end function

end module
