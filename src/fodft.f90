!****h* FHI-aims/fodft
!  NAME
!    fodft
!  SYNOPSIS
module fodft
!  PURPOSE
!  module to handle steps for fragment orbital dft calculations
!  AUTHOR
!    Christoph Schober
!  HISTORY
!    2013
!  USES
    use localorb_io, only: localorb_info
    use timing, only: get_timestamps, get_times, output_times, &
            tot_time_fodft, tot_clock_time_fodft, &
            time_fodft_combine, clock_time_fodft_combine, &
            time_fodft_select, clock_time_fodft_select, &
            time_fodft_lowdin, clock_time_fodft_lowdin, &
            time_fodft_output, clock_time_fodft_output
  implicit none

private

public :: fodft_combine
public :: fodft_select_hab
public :: fodft_set_states
public :: fodft_in_potential
public :: fodft_out_potential

!****** 
contains

!------------------------------------------------------------------------------
!****s* fodft/fodft_combine
!  NAME
!    fodft_combine
!  SYNOPSIS
  subroutine fodft_combine()
!  PURPOSE
!    combines two restart files to a single restart file
!  USES
    use runtime_choices, only: fo_folder1, fo_folder2, occ_fo, fo_type 
    use mpi_tasks, only: myid, mpi_comm_global
    !use physics, only: hamiltonian

    implicit none
!  AUTHOR
!    Christoph Schober
!  HISTORY
!    First version (2015) 
!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE
    character*150 :: info_str

real*8,     dimension(:,:,:), allocatable   :: KS_eigenvalue_work, KS_eigenvalue1, KS_eigenvalue2
real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector1, KS_eigenvector2
real*8,     dimension(:,:,:), allocatable   :: occ_numbers1, occ_numbers2
real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector_sorted
real*8,     dimension(:,:,:), allocatable :: occ_numbers_sorted, occ_numbers
real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector_unsorted
real*8,     dimension(:), allocatable :: occ_work, occ_work_inv
real*8 :: sum_occ

integer, dimension(:), allocatable :: old2new1, new2old1, old2new2, new2old2

integer :: n_basis1, n_states1, n_spin1, i_basis1, i_states1, i_spin1
integer :: n_basis2, n_states2, n_spin2, i_basis2, i_states2, i_spin2
integer :: i_states, i_spin, i_basis, n_basis_work, n_states_work
integer :: n_states_file1, n_states_file2
integer :: n_occ1, n_occ2, n_unocc1, n_unocc2, sum_occ_c, mpierr

call get_timestamps(tot_time_fodft, tot_clock_time_fodft)
call get_timestamps(time_fodft_combine, clock_time_fodft_combine)

! important: We only need all of this on the main processor!
if (myid.eq.0) then 

    ! get everything for fragment 1

    open(file = "../"//fo_folder1//"/restart.frag", unit = 20, status = 'old', form = 'unformatted')
    open(file = "../"//fo_folder1//"/info.frag", unit = 25, status = 'old', form = 'formatted')

    read(25,*) n_occ1
    read(25,*) n_unocc1
    close(unit=25)

    read(20) n_basis1
    read(20) n_states_file1
    read(20) n_spin1

    n_states1 = n_occ1 + n_unocc1

    allocate(KS_eigenvector1(n_basis1,n_states1,n_spin1,1))
    allocate(KS_eigenvalue1(n_states1,n_spin1,1))
    allocate(occ_numbers1(n_states1,n_spin1,1))

    do i_basis1 = 1, n_basis1
        do i_states1 = 1, n_states1
            do i_spin1 = 1, n_spin1
                if (i_states1 <= n_states_file1) then
                    read(20) KS_eigenvector1(i_basis1,i_states1,i_spin1,1)
                else
                    KS_eigenvector1(i_basis1,i_states1,i_spin1,1) = 0.d0
                endif
            end do
        end do
    end do

    do i_states1 = 1, n_states1
        do i_spin1 = 1, n_spin1
            if (i_states1 <= n_states_file1) then
                read(20) KS_eigenvalue1(i_states1,i_spin1,1), occ_numbers1(i_states1,i_spin1,1)
            else
                KS_eigenvalue1(i_states1,i_spin1,1) = 0.d0
                occ_numbers1(i_states1,i_spin1,1) = 0.d0
            endif
        end do
    end do

    close(unit = 20)

    ! get everything for fragment 2
    open(file = "../"//fo_folder2//"/restart.frag", unit = 20, status = 'old', form = 'unformatted')
    open(file = "../"//fo_folder2//"/info.frag", unit = 25, status = 'old', form = 'formatted')

    read(25,*) n_occ2
    read(25,*) n_unocc2
    close(unit=25)

    read(20) n_basis2
    read(20) n_states_file2
    read(20) n_spin2

    n_states2 = n_occ2 + n_unocc2

    allocate(KS_eigenvector2(n_basis2,n_states2,n_spin2,1))
    allocate(KS_eigenvalue2(n_states2,n_spin2,1))
    allocate(occ_numbers2(n_states2,n_spin2,1))

    do i_basis2 = 1, n_basis2
        do i_states2 = 1, n_states2
            do i_spin2 = 1, n_spin2
                if (i_states2 <= n_states_file2) then
                    read(20) KS_eigenvector2(i_basis2,i_states2,i_spin2,1)
                else
                    KS_eigenvector2(i_basis2,i_states2,i_spin2,1) = 0.d0
                endif
            end do
        end do
    end do

    do i_states2 = 1, n_states2
        do i_spin2 = 1, n_spin2
            if (i_states2 <= n_states_file2) then
                read(20) KS_eigenvalue2(i_states2,i_spin2,1), occ_numbers2(i_states2,i_spin2,1)
            else
                KS_eigenvalue2(i_states2,i_spin2,1) = 0.d0
                occ_numbers2(i_states2,i_spin2,1) = 0.d0
            endif
        end do
    end do

    close(unit = 20)

    ! ###### FO-DFT with occupation reset method (similar to CPMD implementation).
    ! if fo_dft from occupied states, reset occupation of fragment 2 spinchannel 2 here
    if ( occ_fo ) then
        allocate(occ_work(n_states2))

        ! For hole-transport we need spinchannel 2, for elec-transport spinchannel 1!
        if (fo_type.eq.'hole'.or.fo_type.eq.'dn') then
            occ_work = reshape(occ_numbers2(1:n_states2, 2:2, 1:1), (/n_states2/))
        elseif (fo_type.eq.'elec'.or.fo_type.eq.'up') then
            occ_work = reshape(occ_numbers2(1:n_states2, 1:1, 1:1), (/n_states2/))
        end if

        allocate(occ_work_inv(n_states2))
        do i_states2 = 1, n_states2
            occ_work_inv(n_states2+1-i_states2) = occ_work(i_states2) 
        end do

        sum_occ = 0.d0
        sum_occ_c = 0
        do i_states2 = 1, n_states2
            if(occ_work_inv(i_states2).NE.0d0) then
                sum_occ = sum_occ + occ_work_inv(i_states2)
                sum_occ_c = sum_occ_c + 1
                occ_work_inv(i_states2) = 0.d0
                if ( 1.d0-sum_occ.LE.1.d-1 ) then
                    write(info_str,'(2X,A)') " "
                    call localorb_info ( info_str )
                    write(info_str,'(2X,A)') "FO-DFT with occupation reset: Resetting HOMO of Fragment2, "
                    call localorb_info ( info_str )
                    write(info_str,'(2X,A, F5.2, A, I1, A)') "removed", sum_occ, " electron from ", sum_occ_c, " state(s)."
                    call localorb_info ( info_str )
                    exit
                end if
            end if
        end do

        do i_states2 = 1, n_states2
            occ_work(i_states2) = occ_work_inv(n_states2+1-i_states2)                                                   
        end do

        if (fo_type.eq.'hole'.or.fo_type.eq.'dn') then
            do i_states2 = 1, n_states2
                occ_numbers2(i_states2, 2, 1) = occ_work(i_states2)
            end do
        elseif (fo_type.eq.'elec'.or.fo_type.eq.'up') then
            do i_states2 = 1, n_states2
                occ_numbers2(i_states2, 1, 1) = occ_work(i_states2)
            end do
        end if
     
        deallocate(occ_work)
        deallocate(occ_work_inv) 
    end if

    ! combine the matrices

    n_basis_work = n_basis1 + n_basis2
    n_states_work = n_states1 + n_states2

    allocate(KS_eigenvector_unsorted(n_basis_work,n_states_work,n_spin1,1))
    KS_eigenvector_unsorted = 0.d0

    allocate(KS_eigenvalue_work(n_states_work,n_spin1,1))
    KS_eigenvalue_work = 0.d0

    allocate(occ_numbers(n_states_work,n_spin1,1))

    ! combine KS_eigenvalues
    KS_eigenvalue_work(1:n_states1, 1:n_spin1, 1:1) = KS_eigenvalue1
    KS_eigenvalue_work(1+n_states1:n_states_work, 1:n_spin1, 1:1) = KS_eigenvalue2

    ! combine occ_numbers
    occ_numbers(1:n_states1, 1:n_spin1, 1:1) = occ_numbers1
    occ_numbers(1+n_states1:n_states_work, 1:n_spin1, 1:1) = occ_numbers2

    ! Sorting is done column-wise
    ! Also sorting for one column only (spin = 1) is irrelevant, since fo-dft always uses 
    ! collinear spins
    ! inserstionsort gives back index array (old2new), use that to sort occupation
    ! numbers and KS_eigenvectors!

    allocate(old2new1(n_states_work))
    allocate(old2new2(n_states_work))
    allocate(new2old1(n_states_work))
    allocate(new2old2(n_states_work))

    call insertionsort(KS_eigenvalue_work(1:n_states_work, 1:1, 1:1), n_states_work, new2old1, old2new1)
    call insertionsort(KS_eigenvalue_work(1:n_states_work, 2:2, 1:1), n_states_work, new2old2, old2new2)

    ! sorting occ_numbers
    allocate(occ_numbers_sorted(n_states_work, n_spin1, 1))

    do i_states = 1, n_states_work
        ! sorting occupation_numbers for spin 1
        occ_numbers_sorted(i_states:i_states, 1:1, 1:1) = &
       &occ_numbers(old2new1(i_states):old2new1(i_states), 1:1, 1:1)

        ! sorting occupation_numbers for spin 2
        occ_numbers_sorted(i_states:i_states, 2:2, 1:1) = &
       &occ_numbers(old2new2(i_states):old2new2(i_states), 2:2, 1:1) 
    end do

    ! now do the combination for KS_eigenvector.. add two m x m matrices
    KS_eigenvector_unsorted(1:n_basis1, 1:n_states1, 1:n_spin1, 1:1) = KS_eigenvector1
    KS_eigenvector_unsorted(1+n_basis1:n_basis_work, 1+n_states1:n_states_work, 1:n_spin1, 1:1) = KS_eigenvector2

    ! write out un-sorted KS_eigenvector for final fo_dft step
    open(file = "fo_ks_eigenvector", unit = 40, status = 'unknown', form = 'unformatted')
    write(40) KS_eigenvector_unsorted
    close(40)

    ! need to sort KS_eigenvector according to KS_eigenvalues!
    allocate(KS_eigenvector_sorted(n_basis_work,n_states_work,n_spin1,1))

    ! get the spin columns of each

    do i_states = 1, n_states_work
        KS_eigenvector_sorted(1:n_basis_work, i_states:i_states, 1:1, 1:1) = &
        & KS_eigenvector_unsorted(1:n_basis_work, old2new1(i_states):old2new1(i_states), 1:1, 1:1)

        KS_eigenvector_sorted(1:n_basis_work, i_states:i_states, 2:2, 1:1) = &
        & KS_eigenvector_unsorted(1:n_basis_work, old2new2(i_states):old2new2(i_states), 2:2, 1:1)
    end do

    ! now write everything to a new restart file
    open(file = "restart.combined", unit = 40, status = 'unknown', form = 'unformatted')

    !first line
    write(40) n_basis_work
    write(40) n_states_work
    write(40) n_spin1  ! n_spin1 == n_spin2 

    ! Write KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
    do i_basis = 1, n_basis_work
        do i_states = 1, n_states_work
            do i_spin = 1, n_spin1
                write(40) KS_eigenvector_sorted(i_basis,i_states,i_spin,1)
            end do
        end do
    end do

    ! Write eigenvalues and occupation numbers.
    do i_states = 1, n_states_work
        do i_spin = 1, n_spin1
            write(40) KS_eigenvalue_work(i_states,i_spin,1), &
                          occ_numbers_sorted(i_states,i_spin,1)
        end do
    end do

    ! Close the file.
    close(unit = 40)

end if ! myid.eq.0

call get_times(tot_time_fodft, tot_clock_time_fodft)
call get_times(time_fodft_combine, clock_time_fodft_combine)

call mpi_barrier(mpi_comm_global, mpierr) 

end subroutine fodft_combine

!------------------------------------------------------------------------------
!****s* fodft/fodft_select_hab
!  NAME
!    fodft_select_hab
!  SYNOPSIS
  subroutine fodft_select_hab(hamiltonian, overlap_matrix)
!  PURPOSE
!  Takes the hamiltonian from scf_solver.f90 and extracts the matrix elements for
!  the state given in the control.in with the fo_orbitals keyword
!  USES
    use runtime_choices, only: fo_folder1, fo_folder2, use_scalapack, use_elpa, &
            fo_type, fo_range1, fo_range2, fo_orb1, fo_orb2, fo_verbosity
    use scalapack_wrapper, only: mxld, mxcol, ham, ovlp, &
            construct_hamiltonian_scalapack, set_full_matrix_real, &
            get_scalapack_global_rmatrix
    use mpi_utilities, only: myid
    use dimensions, only: n_basis, n_states, n_spin
    use constants, only: Hartree
    use linalg, only: lowdin_lapack, lowdin_elpa, get_rinverse_elpa, &
            get_rinverse_lapack, simple_matmul
    use aims_memory_tracking, only : aims_deallocate
    implicit none
!  AUTHOR
!    Christoph Schober
!  INPUTS
!    o hamiltonian
!  OUTPUT
!    none
!  SEE ALSO
!    FO-DFT publication. 
!  SOURCE
real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector_unsorted
real*8,     dimension(:,:), allocatable :: overlap_work, hamiltonian_work, KS_ev_work, transformed_mat, X_mat, KS_ev_ortho, ovlp_trans, ham_trans, hamiltonian, tmp
real*8,      dimension(:), allocatable :: overlap_matrix

real*8 ::  v_ab, v_aa, v_bb, s_ab
real*8, dimension(:,:), allocatable :: h_ab, h_ab_full

integer :: i_basis, k, i_range, j_range, spin_channel
integer :: n_states_work, n_basis1, n_states1, n_spin1
integer :: n_basis2, n_states2, n_spin2

character*150 :: info_str

call get_timestamps(tot_time_fodft, tot_clock_time_fodft)
call get_timestamps(time_fodft_select, clock_time_fodft_select)

write(info_str,'(A)') " "
call localorb_info ( info_str )
write(info_str,'(A)') "  --------------------------------------------------------"
call localorb_info ( info_str )
write(info_str,'(A)') "  Calculation of electronic coupling values with FO-DFT"
call localorb_info ( info_str )
write(info_str,'(A)') "  --------------------------------------------------------"
call localorb_info ( info_str )

! get number of basis, states and the unsorted eigenvector for both fragments
open(file = "../"//fo_folder1//"/restart.frag", unit = 25, status = 'old', form = 'unformatted')
open(file = "../"//fo_folder2//"/restart.frag", unit = 35, status = 'old', form = 'unformatted')

read(25) n_basis1
read(25) n_states1
read(25) n_spin1

read(35) n_basis2
read(35) n_states2
read(35) n_spin2

n_states_work = n_states
close(unit = 25)
close(unit = 35)

open(file = "fo_ks_eigenvector", unit = 40, status = 'unknown', form = 'unformatted')
allocate(KS_eigenvector_unsorted(n_basis, n_states_work, n_spin, 1))
read(40) KS_eigenvector_unsorted
close(40)
if (fo_type.eq.'hole'.or.fo_type.eq.'dn') then
    spin_channel = 2
elseif (fo_type.eq.'elec'.or.fo_type.eq.'up') then
    spin_channel = 1
end if


allocate(KS_ev_work(n_basis, n_states_work))

if (fo_type.eq.'hole'.or.fo_type.eq.'dn') then
    KS_ev_work = reshape(KS_eigenvector_unsorted(1:n_basis, 1:n_states_work, 2:2, 1:1), (/n_basis, n_states_work/))
elseif (fo_type.eq.'elec'.or.fo_type.eq.'up') then
    KS_ev_work = reshape(KS_eigenvector_unsorted(1:n_basis, 1:n_states_work, 1:1, 1:1), (/n_basis, n_states_work/))
end if

deallocate(KS_eigenvector_unsorted)

!if (packed_matrix_format /= PM_none) then
if (use_scalapack.and.use_elpa) then
    call construct_hamiltonian_scalapack( hamiltonian )
    allocate(hamiltonian_work(mxld, mxcol))
    hamiltonian_work = ham(:,:,spin_channel)
    call aims_deallocate(ham, "ham") !!!! This is (was) the global AIMS hamiltonian.
    call set_full_matrix_real(hamiltonian_work)

else ! = LAPACK
    allocate(hamiltonian_work(n_basis, n_basis))
    if (fo_type.eq.'hole'.or.fo_type.eq.'dn') then
        do i_basis = 1, n_basis
            k = i_basis*(i_basis-1)/2
            hamiltonian_work(1:i_basis, i_basis) = hamiltonian(1+k:i_basis+k, 2)
            hamiltonian_work(i_basis, 1:i_basis) = hamiltonian(1+k:i_basis+k, 2)
        end do
    elseif (fo_type.eq.'elec'.or.fo_type.eq.'up') then
        do i_basis = 1, n_basis
            k = i_basis*(i_basis-1)/2
            hamiltonian_work(1:i_basis, i_basis) = hamiltonian(1+k:i_basis+k, 1)
            hamiltonian_work(i_basis, 1:i_basis) = hamiltonian(1+k:i_basis+k, 1)
        end do
    end if ! hole vs elec
end if

! At this point we can distinguish between two paths: LAPACK and scaLAPACK.
!
! For LAPACk we have now:
! o KS_ev_work(n_basis, n_states)
! o hamiltonian_work(n_basis, n_basis)
! o overlap_matrix(n_basis*(n_basis+1)/2)
!
! All subsequent calculations are done with LAPACK routines for m*n (dgemm, dpotrf, dpotri)
! or n*(n+1)/2 (for diagonalize_overlap).
!
! For scaLAPACK we have:
! o KS_ev_work(mxld, mxcol)  [CHECK IF THIS IS CORRECT FOR m*n if m != n]
! o hamiltonian_work(mxld, mxcol)
! o overlap_work(mxld, mxcol)
!
! All subsequent calculations are done with scaLAPACK or ELPA routines for
! 2D-block cyclic distribution (pdgemm, invert_trm_real, etc)
!
!PSEUDOCODE
!
! LOWDIN
! U**T*(s**-0.5)*U = S**-0.5 = X
! 1) s, U = diag(S)
! 2) s**-0.5 = diag_mat(s**-0.5)
! 3) us = dgemm(U, s**-0.5)
! 4) X = S**-0.5 = dgemm(us, U)
call get_timestamps(time_fodft_lowdin, clock_time_fodft_lowdin)

if (use_scalapack) then
    allocate(X_mat(mxld, mxcol))
    X_mat = 0.d0
    call lowdin_elpa(ovlp, X_mat, n_basis)
else
    allocate(overlap_work(n_basis, n_basis))
        do i_basis = 1, n_basis
            k = i_basis*(i_basis-1)/2
            overlap_work(1:i_basis, i_basis) = overlap_matrix(1+k:i_basis+k)
            overlap_work(i_basis, 1:i_basis) = overlap_matrix(1+k:i_basis+k)
        end do

    allocate(X_mat(n_basis, n_basis))
    X_mat = 0.d0
    !call linalg_lowdin(overlap_matrix, X_mat)
    call lowdin_lapack(overlap_work, X_mat)
end if

call get_times(time_fodft_lowdin, clock_time_fodft_lowdin)

! H' = X**T*H*X
! 1) xh = dgemm(X, H)
! 2) H' = dgemm(xh, X)
allocate(tmp(size(hamiltonian_work, 1), size(hamiltonian_work, 2)))
if (use_scalapack) then
    call simple_matmul('T', 'N', X_mat, hamiltonian_work, tmp, n_basis, n_basis, n_basis)
    call simple_matmul('N', 'N', tmp, X_mat, hamiltonian_work, n_basis, n_basis, n_basis)
else
    call dgemm('T', 'N', n_basis, n_basis, n_basis, 1.d0, X_mat, ubound(X_mat, 1), hamiltonian_work, ubound(hamiltonian_work, 1), 0.d0, tmp, ubound(tmp, 1))
    call dgemm('N', 'N', n_basis, n_basis, n_basis, 1.d0, tmp, ubound(tmp, 1), X_mat, ubound(X_mat, 1), 0.d0, hamiltonian_work, ubound(hamiltonian_work, 1))
end if

deallocate(tmp)

! X**-1 = invert(X)
!allocate(inverse_X(size(hamiltonian_work, 1), size(hamiltonian_work, 2)))
!inverse_X = X_mat

if (use_scalapack) then
    call get_rinverse_elpa(X_mat, n_basis)
else
    call get_rinverse_lapack(X_mat, n_basis)
end if

! Continue with LAPACK matrices to avoid the full scalapack setup for arrays different from n_basis x n_basis
! To do this, X_mat, hamiltonian_work, ... need to be collected from distributed storage...
!
! If memory in the dimer step get a problem, this should be changed to
! support full scalapack arrays.
if (use_elpa) then
    allocate(transformed_mat(n_basis, n_basis))
    call get_scalapack_global_rmatrix( X_mat, transformed_mat)
    deallocate(X_mat)
    allocate(X_mat(n_basis, n_basis))
    X_mat = transformed_mat
    deallocate(transformed_mat)

    allocate(transformed_mat(n_basis, n_basis))
    call get_scalapack_global_rmatrix( hamiltonian_work, transformed_mat)
    deallocate(hamiltonian_work)
    allocate(hamiltonian_work(n_basis, n_basis))
    hamiltonian_work = transformed_mat
    deallocate(transformed_mat)

end if

! C' (C_ortho) = dgemm(X**-1, C)

allocate(KS_ev_ortho(size(KS_ev_work, 1), size(KS_ev_work, 2)))

call simple_matmul('N', 'N', X_mat, KS_ev_work, KS_ev_ortho)

deallocate(KS_ev_work)

! H_trans = C'**T*H'*C'
! 1) ch = dgemm(C'**T, H)
! 2) H_trans = dgemm(ch, C')

allocate(tmp(size(KS_ev_ortho, 2), size(KS_ev_ortho, 1)))
call simple_matmul('T', 'N', KS_ev_ortho, hamiltonian_work, tmp)
deallocate(hamiltonian_work)

allocate(ham_trans(size(KS_ev_ortho, 2), size(KS_ev_ortho, 2)))
call simple_matmul('N', 'N', tmp, KS_ev_ortho, ham_trans)
deallocate(tmp)


! S_trans = dgemm(C'**T, C')
allocate(ovlp_trans(size(KS_ev_ortho, 2), size(KS_ev_ortho, 2)))
call simple_matmul('T', 'N', KS_ev_ortho, KS_ev_ortho, ovlp_trans)
deallocate(KS_ev_ortho)

! From here on only ovlp_trans and ham_trans are needed!

allocate(h_ab(fo_range1, fo_range2))

write(info_str,'(2X,A)') "Collecting information for output..."
call localorb_info ( info_str )

call get_timestamps(time_fodft_output, clock_time_fodft_output)

allocate(h_ab_full(n_states1, n_states2))

do i_range = 1, n_states1
    do j_range = 1, n_states2
        v_aa = ham_trans(i_range, i_range)
        v_bb = ham_trans(n_states1+j_range, n_states1+j_range)

        v_ab = ham_trans(i_range, n_states1+j_range)
        s_ab = ovlp_trans(i_range, n_states1+j_range)

        ! final V_ab
        h_ab_full(i_range, j_range) = (1-s_ab**(2))**(-1)*(v_ab-s_ab*(v_aa+v_bb)/2)
    end do
end do

! Write out full 
if (fo_verbosity.gt.0) then
    if (myid.eq.0) then
        open(file = "full_hab_submatrix", unit = 20, status = 'unknown', form = 'formatted')
        write(20,*) h_ab_full
        close(20)
    end if
end if

do i_range = 1, fo_range1
    do j_range = 1, fo_range2
        ! calculate v_ab_upper AND v_ab_lower ... triangle
        v_aa = ham_trans(fo_orb1+i_range-1, fo_orb1+i_range-1)
        v_bb = ham_trans(n_states1+fo_orb2+j_range-1, n_states1+fo_orb2+j_range-1)

        v_ab = ham_trans(fo_orb1+i_range-1, n_states1+fo_orb2+j_range-1)
        s_ab = ovlp_trans(fo_orb1+i_range-1, n_states1+fo_orb2+j_range-1)

        ! final V_ab
        h_ab(i_range, j_range) = (1-s_ab**(2))**(-1)*(v_ab-s_ab*(v_aa+v_bb)/2)

        if (fo_verbosity.eq.2) then
            write(info_str,'(2X,A)') "Verbose output selected. Writing values for v_ab, s_ab, v_aa and v_bb"
            call localorb_info ( info_str )

            write(info_str,'(2X,A)') " "
            call localorb_info ( info_str )
            write(info_str,'(2X,A, I4, 2X,A, I4)') "States:", fo_orb1+i_range-1, "->", fo_orb2+j_range-1
            call localorb_info ( info_str )
            write(info_str,'(2X,A, F10.8, 2X,A)') "v_aa: ", v_aa, "Ha"
            call localorb_info ( info_str )
            write(info_str,'(2X,A, F10.8, 2X,A)') "v_bb: ", v_bb, "Ha"
            call localorb_info ( info_str )
            write(info_str,'(2X,A, F10.8, 2X,A)') "s_ab: ", s_ab, "Ha"
            call localorb_info ( info_str )
            write(info_str,'(2X,A, F10.8, 2X,A)') "v_ab: ", v_ab, "Ha"
            call localorb_info ( info_str )
        end if !fo_verbosity

    end do
end do

write(info_str,'(2X,A, A, A)') "H_ab for ", fo_type, "-transfer"
call localorb_info ( info_str )
write(info_str,'(2X,A, T35,I3, T40,A, T44,I3)') "States evaluated for Fragment 1:", fo_orb1, "to", fo_orb1+fo_range1-1
call localorb_info ( info_str )
write(info_str,'(2X,A, T35,I3, T40,A, T44,I3)') "States evaluated for Fragment 2:", fo_orb2, "to", fo_orb2+fo_range2-1
call localorb_info ( info_str )
write(info_str,'(2X,A)') "Please remember, these states are numbered according to the fragments,"
call localorb_info ( info_str )
write(info_str,'(2X,A)') "NOT the combined system!"
call localorb_info ( info_str )
!call localorb_info ( info_str )
write(info_str,'(1X,A)') " "
call localorb_info ( info_str )


write(info_str,'(2X,A, T19,A)') "st1 -> st2:", "h_ab"
call localorb_info ( info_str )
write(info_str,'(1X,A)') " -----------  ----------   ---"
call localorb_info ( info_str )
do i_range = 1, fo_range1
    do j_range = 1, fo_range2
        write(info_str,'(T3,I3.3, T7,A, T10,I3.3, T13,A, T16,F10.2, T29,A)')&
                                 fo_orb1+i_range-1, "->", fo_orb2+j_range-1, &
                                 ":", h_ab(i_range, j_range)*Hartree*1000, "meV"
        call localorb_info ( info_str )
        write(info_str,'(2X,A, T13,A, T16,F10.3, T29,A)') ".", ":", h_ab(i_range, j_range)*1000, "mH"
        call localorb_info ( info_str )
    end do
end do
call get_times(time_fodft_output, clock_time_fodft_output)

call get_times(tot_time_fodft, tot_clock_time_fodft)
call get_times(time_fodft_select, clock_time_fodft_select)

write(info_str,'(2X,A)') ""
call localorb_info ( info_str )

call output_times("2X", "Total time to calculate H_ab", time_fodft_select, clock_time_fodft_select)
write(info_str,'(2X,A)') "-------------------------------------------------------------------------------"
call localorb_info ( info_str )
call output_times("2X", "Time for Loewdin orthogonalisation", time_fodft_lowdin, clock_time_fodft_lowdin)
call output_times("2X", "Time for H_ab submatrix evaluation", time_fodft_output, clock_time_fodft_output)


write(info_str,'(2X,A)') "-------------------------------------------------------------------------------"
call localorb_info ( info_str )

write(info_str,'(2X,A)') "Leaving FO-DFT routine."
call localorb_info ( info_str )


end subroutine fodft_select_hab

subroutine fodft_set_states(n_empty)
!  PURPOSE
! Get the number of empty states for the fragment calculations
!  USES
  use mpi_tasks, only: myid, mpi_comm_global, use_mpi, aims_stop
  use mpi_utilities, only: MPI_INTEGER
  use runtime_choices, only: fo_folder1, fo_folder2
  
  implicit none

  integer :: n_empty, mpierr
  integer :: empty_states1, empty_states2
  logical :: file1_exists, file2_exists

if (myid.eq.0) then
    ! this assumes that in 'min(n_states, n_states_occupied + 3)' for each fragment
    ! the n_states_occupied + 3 is smallest. 
    ! Will break for very small systems.

    inquire(FILE="../"//fo_folder1//"/info.frag", EXIST=file1_exists)
    inquire(FILE="../"//fo_folder2//"/info.frag", EXIST=file2_exists)

    if (file1_exists) then
        open(file = "../"//fo_folder1//"/info.frag", unit = 25, status = 'old', form = 'formatted')
        read(25,*) empty_states1
        read(25,*) empty_states1
        close(unit = 25)
    else
        call aims_stop("No 'info.frag' file for FODFT fragment 1 found. Please add file manually (see docs) &
             &or use option 'fodft fragment' in fragment calculation")
    endif

    if (file2_exists) then
        open(file = "../"//fo_folder2//"/info.frag", unit = 25, status = 'old', form = 'formatted')
        read(25,*) empty_states2
        read(25,*) empty_states2
        close(unit = 25)
    else
        call aims_stop("No 'info.frag' file for FODFT fragment 2 found. Please add file manually (see docs) &
             &or use option 'fodft fragment' in fragment calculation")
    endif

    n_empty = empty_states1 + empty_states2
end if ! myid.eq.0

if (use_mpi) then
    call mpi_bcast(n_empty, 1, MPI_INTEGER, 0, mpi_comm_global, mpierr )
    call mpi_barrier(mpi_comm_global, mpierr)
end if

end subroutine fodft_set_states

subroutine fodft_out_potential( local_potential_parts, output_mode, local_fo_potential ) 
!
! subroutine to write out the (hartree)potential on the native aims grid for all points and mpi-processes
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use mpi_tasks
  use mpi_utilities
  use species_data

  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: local_potential_parts
!  local variables

  real*8  :: grid_coord(3)

!  counters

  integer              :: i_coord
  integer              :: output_mode
  integer              :: i_point, i_my_batch, i_index
  character( len = 8 ) :: myid_string
  character*150 :: info_str
  real*8, optional     :: local_fo_potential(n_full_points)

  logical         :: file_exists

  call get_my_task()
  
  ! convert mpi task id to a string
  write(unit=myid_string, fmt='(i8.8)') myid 

  inquire(FILE="output_grid_"//myid_string, EXIST=file_exists)
  
  if (.not.file_exists) then
    open (50, file="output_grid_"//myid_string, form='unformatted', status='new')
  
    write(50) n_full_points

    i_point = 0
    ! loop over all batches
    do i_my_batch = 1, n_my_batches, 1
  
      ! loop over all points in a batch
      do i_index = 1, batches(i_my_batch)%size, 1         
          i_point = i_point + 1    
          grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)*bohr
          write(50) (grid_coord(i_coord), i_coord=1,3,1), local_potential_parts(i_point)!, partition_tab_std(i_point)
      end do
     
    end do
    close(unit=50)
  else
    open (50, file="output_grid_"//myid_string//"_new", form='unformatted', status='new')
    write(50) n_full_points

    ! subtract the additional (from embedded fragment) hartree potential    
    local_potential_parts = local_potential_parts - local_fo_potential

    i_point = 0
    ! loop over all batches
    do i_my_batch = 1, n_my_batches, 1
  
      ! loop over all points in a batch
      do i_index = 1, batches(i_my_batch)%size, 1
          i_point = i_point + 1
          grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)*bohr
          write(50) (grid_coord(i_coord), i_coord=1,3,1), local_potential_parts(i_point)!, partition_tab_std(i_point)
      end do

    end do
    close(unit=50)

  end if !(.not.file_exists)
  
  write(info_str,'(A)') "  ****************** FO_DFT embedding ************************"
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "FO-DFT fragment calculation with embedding potential."
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "Writing Hartree potential on AIMS-grid to file 'output_grid_MPI_TASK'.."
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "Use this file for embedding of the potential into the other FO-DFT fragment."
  call localorb_info ( info_str )
  write(info_str,'(A)') "  ************************************************************"
  call localorb_info ( info_str )

end subroutine fodft_out_potential

subroutine fodft_in_potential( local_fo_potential )

! Routine to read in saved hartree potential on aims grid
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use mpi_tasks
  use mpi_utilities
  use species_data

implicit none

!  ARGUMENTS

  character( len = 8 ) :: myid_string
  real*8  :: grid_coord(3), fo_coord(3)
  real*8, dimension(:), allocatable  :: fo_potential, &
                                     local_fo_potential, fo_pot_tosync, &
                                     fo_pot_synced, fo_coord_synced
  real*8, dimension(:,:), allocatable :: fo_grid_coord, fo_coord_tosync, new_fo_coord
!  real*8, dimension(:,:), allocatable  :: fo_grid_coord, fo_potential
!  counters

  integer              :: i_coord, n_fo_points, mask_counter, loop_counter, source
  integer              :: i_pot, i_my_batch, i_index, mpierr, dest, &
              new_temp_fo_points, mpi_status(MPI_STATUS_SIZE), &
                sendtag, recvtag, total_equal, total_unequal, total_points
  integer              :: i_point, equal_counter, unequal_counter, temp_fo_points
  logical              :: all_assigned, recme, imdone, alldone(n_tasks)
  !integer              :: sender_list(n_tasks), receiver_list(n_tasks)

  character*150 :: info_str
  call get_my_task()
  
  ! convert mpi task id to a string
  write(unit=myid_string, fmt='(i8.8)') myid 
  
  open(50, file="output_grid_"//myid_string, form='unformatted', status='old')
  read(50) n_fo_points

  if (.not.allocated(fo_potential))     allocate(fo_potential(n_fo_points))
  if (.not.allocated(fo_grid_coord)) allocate(fo_grid_coord(3, n_fo_points))
  do i_point = 1, n_fo_points, 1
      read(50) (fo_grid_coord(i_coord, i_point), i_coord=1,3,1), fo_potential(i_point)  
  end do
  
  close(unit=50)

  write(info_str,'(A)') "  ****************** FO_DFT embedding ************************"
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "FO-DFT fragment calculation with embedding potential."
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "Found Hartree-potential from other fragment."
  call localorb_info ( info_str )
  write(info_str,'(2X,A)') "Started to re-assign grid points for external embedding potential..."
  call localorb_info ( info_str )

  all_assigned = .false.

  !i_point = 0
  equal_counter = 0
  unequal_counter = 0
  loop_counter = 0
  recme = .true.
  do while (.not.all_assigned)
  ! loop until all loaded potential points are assigned to some MPI task..
 
  i_point = 0
  loop_counter = loop_counter + 1
  mask_counter = 0

      if (loop_counter.gt.1) then !re-adjust dimensions for fo_*_tosync, only applies for MPI runs with ntasks > 1

        n_fo_points = new_temp_fo_points 
         
        if (allocated(fo_grid_coord))   deallocate(fo_grid_coord)
        if (allocated(fo_potential))    deallocate(fo_potential)
        
        allocate(fo_grid_coord(3, n_fo_points))
        allocate(fo_potential(n_fo_points))

        fo_grid_coord = new_fo_coord
        fo_potential = fo_pot_synced

      end if !(loop_counter)
      ! loop over fo_potential
  
    ! loop over all batches
    do i_my_batch = 1, n_my_batches, 1
  
      ! loop over all points in a batch
      do i_index = 1, batches(i_my_batch)%size, 1         
      i_point = i_point + 1    
      grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)*bohr

       do i_pot = 1, n_fo_points, 1
              fo_coord(:) = fo_grid_coord(:, i_pot)
 
         if ((abs(grid_coord(1)-fo_coord(1))<0.00001).and.&
             (abs(grid_coord(2)-fo_coord(2))<0.00001).and.&
             (abs(grid_coord(3)-fo_coord(3))<0.00001)) then
     
            equal_counter = equal_counter + 1
            mask_counter = mask_counter + 1

            do i_coord = 1, 3, 1
                fo_grid_coord(i_coord, i_pot) = 123456789.0
            end do

            local_fo_potential(i_point) = fo_potential(i_pot)
            fo_potential(i_pot) = 123456789.0
          else
            unequal_counter = unequal_counter + 1
          end if
  
      end do ! i_index
    end do !i_batch
  end do ! i_pot

  if (.not.use_mpi) then
    total_points = n_full_points
  end if

  ! Now, do the MPI part:
  !     1. Book-keeping: - Which tasks still have unassigned points in fo_potential
  !                      - Which tasks still have unmatched points in local_fo_potential
  if (use_mpi) then

      temp_fo_points = n_fo_points-mask_counter

      if (allocated(fo_pot_tosync)) deallocate(fo_pot_tosync)
      allocate(fo_pot_tosync(temp_fo_points))
      if (allocated(fo_coord_tosync)) deallocate(fo_coord_tosync)
      allocate(fo_coord_tosync(3, temp_fo_points))

          mask_counter = 0
          do i_pot = 1, n_fo_points
                  if (.not.fo_potential(i_pot).eq.123456789.0) then
                      mask_counter = mask_counter + 1
                      fo_pot_tosync(mask_counter) = fo_potential(i_pot)
                      fo_coord_tosync(:, mask_counter) = fo_grid_coord(:, i_pot)
                  end if
                  
          end do

     ! get send destination depending on myid 
    if (myid+1.lt.n_tasks) then
        dest = myid + 1
    else
        dest = 0
    end if

    if (myid.eq.0) then
        source = n_tasks - 1
    else
        source = myid - 1
    end if


    ! only do the following part when we are in fact using MPI / more than 1 core.

        ! send size of fo_pot_tosync to myid+1
        ! receive size of previous fo_pot_tosync from myid-1
        sendtag = 111
        recvtag = 111
        call mpi_sendrecv(temp_fo_points, 1, MPI_INTEGER, dest, sendtag, &
                          new_temp_fo_points, 1, MPI_INTEGER, source, recvtag, mpi_comm_global, mpi_status, mpierr)

        if (allocated(fo_pot_synced)) deallocate(fo_pot_synced)
        allocate(fo_pot_synced(new_temp_fo_points))
        if (allocated(fo_coord_synced)) deallocate(fo_coord_synced)
        allocate(fo_coord_synced(3*new_temp_fo_points))
    
        call mpi_sendrecv(fo_pot_tosync, temp_fo_points, MPI_REAL8, dest, sendtag, &
                          fo_pot_synced, new_temp_fo_points, MPI_REAL8, source, recvtag, mpi_comm_global, mpi_status, mpierr)
        call mpi_sendrecv(fo_coord_tosync, 3*temp_fo_points, MPI_REAL8, dest, sendtag, &
                          fo_coord_synced, 3*new_temp_fo_points, MPI_REAL8, source, recvtag, mpi_comm_global, mpi_status, mpierr)
     
        !reshape the coords array to correct form
        if (allocated(new_fo_coord)) deallocate(new_fo_coord)
        allocate(new_fo_coord(3, new_temp_fo_points))

        new_fo_coord = reshape(fo_coord_synced, (/3, new_temp_fo_points/))

        ! check if all points on all cores are assigned
        if ((temp_fo_points.eq.0).and.(ubound(local_fo_potential,1).eq.n_full_points)) then
            imdone = .true.
        else
            imdone = .false.
        end if

        call mpi_gather(imdone, 1, MPI_LOGICAL, alldone(myid+1), 1, MPI_LOGICAL, 0, mpi_comm_global, mpierr)

        call mpi_barrier(mpi_comm_global, mpierr)

        if (myid.eq.0) then
            if (ALL(alldone, 1)) then
                all_assigned = .true.
            end if        
        end if

        call mpi_bcast(all_assigned, 1, MPI_LOGICAL, 0, mpi_comm_global, mpierr )
        
        call mpi_barrier(mpi_comm_global, mpierr)

    else !(= serial calculation, if nothing gone wrong all points are assigned in first loop)
        all_assigned = .true.

    end if !(use mpi)


    end do ! while .not.all_assigned
      
    if (use_mpi) then     
          ! now do some statistics
          call mpi_reduce(equal_counter, total_equal, 1, MPI_INTEGER, MPI_SUM, 0, mpi_comm_global, mpierr)
          call mpi_reduce(unequal_counter, total_unequal, 1, MPI_INTEGER, MPI_SUM, 0, mpi_comm_global, mpierr)
          call mpi_reduce(n_full_points, total_points, 1, MPI_INTEGER, MPI_SUM, 0, mpi_comm_global, mpierr)

          ! to make the output valid even for serial calculations: 
          equal_counter = total_equal
          unequal_counter = total_unequal

    end if !(use mpi)

  write(info_str,'(2X,A, I10, 1X,A, I8, 1X,A)') "Successfully assigned ", equal_counter, "of ", total_points, "grid points."
  call localorb_info ( info_str )
  write(info_str,'(2X,A, I20, 1X,A)') "Needed ", equal_counter+unequal_counter, "steps..."
  call localorb_info ( info_str )
  write(info_str,'(A)') "  ************************************************************"
  call localorb_info ( info_str )

! Deallocations for all the used arrays that are not needed anymore

if (allocated(fo_pot_tosync)) deallocate(fo_pot_tosync)
if (allocated(fo_pot_synced))  deallocate(fo_pot_synced)
if (allocated(fo_coord_synced))  deallocate(fo_coord_synced)
if (allocated(fo_grid_coord))   deallocate(fo_grid_coord)
if (allocated(fo_potential))    deallocate(fo_potential)
if (allocated(fo_coord_tosync)) deallocate(fo_coord_tosync)
if (allocated(new_fo_coord)) deallocate(new_fo_coord)

end subroutine fodft_in_potential

end module fodft
