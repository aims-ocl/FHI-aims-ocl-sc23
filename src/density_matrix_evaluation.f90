!****h* FHI-aims/density_matrix_evaluation
!  NAME
!    density_matrix_evaluation
!  SYNOPSIS

module density_matrix_evaluation

  !  PURPOSE
  !    Provide procedures for evaluating the density matrix.
  !  USES
  implicit none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  private

  public :: evaluate_densmat
  public :: evaluate_k_densmat
  public :: accumulate_k_densmat
  public :: output_densmat

contains

  !****s* density_matrix_evaluation/evaluate_densmat
  !  NAME
  !    evaluate_densmat
  !  SYNOPSIS
  subroutine evaluate_densmat &
  ( KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
  density_matrix, density_matrix_sparse, i_spin, force_packed &
  )
    !  PURPOSE
    !    Evaluates the density matrix
    !  USES
    use aims_memory_tracking, only : aims_allocate, aims_deallocate
    use dimensions, only: n_basis, n_states, n_spin, n_k_points_task, &
        n_k_points, n_centers_basis_T, n_hamiltonian_matrix_size
    use elsi_wrapper, only: rwh_r, rwh_w, aims_elsi_write_mat
    use load_balancing, only: batch_perm, use_batch_permutation, &
        get_full_local_matrix
    use localorb_io, only: use_unit, OL_norm, localorb_info
    use mpi_tasks, only: myid, n_tasks, aims_stop
    use pbc_lists, only: k_weights
    use rel_x2c_mod, only: dim_matrix_rel
    use restart_elsi, only: elsi_restart_scalapack, elsi_restart_lapack
    use runtime_choices, only: packed_matrix_format, pm_none, use_scalapack, &
        real_eigenvectors, use_local_index, flag_rel,rel_x2c, rel_4c_dks, &
        elsi_read_dm, elsi_read_dm_done, do_elsi_rw, elsi_write_dm_this, &
        pm_index
    use scalapack_wrapper, only: ham, ham_complex, my_k_point, &
        construct_dm_scalapack, get_sparse_local_matrix_scalapack, &
        get_sparse_matrix_scalapack, get_full_matrix_scalapack, &
        set_full_matrix_real, set_full_matrix_complex
    use synchronize_mpi, only: sync_density_matrix_sparse, sync_density_matrix,&
                               sync_vector_complex
    use timing_core, only: get_times, get_timestamps
    use timing, only: tot_time_matrix_io, tot_clock_time_matrix_io
    implicit none

    !  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector
    complex*16 :: KS_eigenvector_complex(*) ! (n_basis, n_states, n_spin,n_k_points_task)
                                            ! or (2*dim_matrix_rel, n_states, n_spin,n_k_points_task)
    real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
    ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
    !        properly k-weighted (and are thus not the "true" occupation numbers.)
    !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
    integer, intent(IN) :: i_spin
    real*8, dimension(n_centers_basis_T,n_centers_basis_T), intent(OUT) :: density_matrix
    ! when this routine is called, density_matrix_sparse has either the dimension
    ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
    ! so we declare it here as a 1D assumed size array
    real*8, intent(OUT) :: density_matrix_sparse(*)
    logical, intent(IN) :: force_packed

    !  INPUTS
    !   o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !   o occ_numbers -- occupation of states, with k-weighting already applied
    !   o i_spin -- spin index
    !   o force_packed -- use packed storage (Lapack-type in ..._sparse) even if
    !                   packed_matrix_format is PM_none
    !
    !  OUTPUT
    !   o density_matrix -- density matrix if non-packed matrix is in use
    !   o density_matrix_sparse -- density matrix if packed matrix is in use
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !    Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    !

    !     other local variables
    real*8, dimension(:,:), allocatable :: kdm
    complex*16, dimension(:,:), allocatable :: kdm_complex
    complex*16, dimension(:,:,:), allocatable :: kdm_complex_rel
    real*8, dimension(:), allocatable :: density_matrix_temp
    integer :: i_k_point, i_k, n
    integer :: info, i
    character(*), parameter :: func = "evaluate_densmat"

    ! ELSI restart
    character*100 :: mat_file
    real*8 :: n_electrons_file
    integer :: n_basis_file
    integer :: mxld_file
    integer :: mxcol_file

    ! Timing output
    real*8 :: t_read
    real*8 :: t_write
    real*8 :: t_read1
    real*8 :: t_write1
    real*8 :: tmp
    character*200 :: info_str

    call localorb_info("Evaluating density matrix",use_unit,"(2X,A)",OL_norm)

    if(packed_matrix_format /= PM_none)then
       if(use_batch_permutation > 0) then
          n = batch_perm(use_batch_permutation)%n_local_matrix_size
          density_matrix_sparse(1:n) = 0.d0
       else
          density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.d0
       endif
    else if (force_packed) then
       density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.d0
    else
       density_matrix = 0.d0
    end if

    if(use_scalapack)then
       if(elsi_read_dm .and. do_elsi_rw) then
          ! Could have been done in calculate_fock_matrix_p0.f90
          if(elsi_read_dm_done) then
             if(real_eigenvectors) then
                ham(:,:,i_spin) = ham(:,:,i_spin)*k_weights(my_k_point)
             else
                ham_complex(:,:,i_spin) = ham_complex(:,:,i_spin)&
                   *k_weights(my_k_point)
             endif
          else
             call get_timestamps(tmp,t_read)

             if(real_eigenvectors) then
                call elsi_restart_scalapack(i_spin,ham(:,:,i_spin))
             else
                call elsi_restart_scalapack(i_spin,ham_complex(:,:,i_spin))
             endif

             call get_times(tmp,t_read,tot_time_matrix_io,&
                  tot_clock_time_matrix_io)
          endif
       else
          call construct_dm_scalapack(occ_numbers,i_spin)
       endif

       if(elsi_write_dm_this .and. do_elsi_rw) then
          call get_timestamps(tmp,t_write)

          write(mat_file,"(A,I2.2,A,I6.6,A)") &
             "D_spin_",i_spin,"_kpt_",my_k_point,".csc"

          if(real_eigenvectors) then
             call set_full_matrix_real(ham(:,:,i_spin))
             call aims_elsi_write_mat(rwh_w,trim(mat_file),ham(:,:,i_spin))
          else
             call set_full_matrix_complex(ham_complex(:,:,i_spin))
             call aims_elsi_write_mat(rwh_w,trim(mat_file),&
                  ham_complex(:,:,i_spin))
          endif

          call get_times(tmp,t_write,tot_time_matrix_io,&
               tot_clock_time_matrix_io)
       endif

       select case (packed_matrix_format)
       case(PM_index)
          if(use_local_index) then
             if(use_batch_permutation > 0) then
                call get_full_local_matrix(density_matrix_sparse, i_spin)
             else
                call get_sparse_local_matrix_scalapack(density_matrix_sparse, i_spin)
             endif
          else
             call get_sparse_matrix_scalapack(density_matrix_sparse,i_spin)
             call sync_density_matrix_sparse(density_matrix_sparse)
          endif

       case(PM_none)
          if (force_packed) then
             call aims_stop("LAPACK-packed matrix not for ScaLAPACK", func)
          end if
          call get_full_matrix_scalapack( density_matrix, i_spin )
          call sync_density_matrix(density_matrix)
       case default
          call aims_stop("Invalid packing!", func)
       end select

    else

       !------------ construct density matrix -------------------------------------------

       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
          call aims_allocate( kdm_complex_rel, 2*dim_matrix_rel, 2*dim_matrix_rel, &
                              n_k_points, "kdm_complex_rel" )
          kdm_complex_rel = (0.d0,0.d0)
          i_k = 0
          do i_k_point = 1,n_k_points
             if(myid == mod(i_k_point,n_tasks) .and. myid <= n_k_points)then
                i_k = i_k+1
                call evaluate_k_densmat_rel(kdm_complex_rel(1,1,i_k_point),occ_numbers, &
                     KS_eigenvector_complex,i_spin,i_k_point,i_k)
             endif
          enddo
          call sync_vector_complex(kdm_complex_rel,2*dim_matrix_rel*2*dim_matrix_rel*n_k_points)
          call accumulate_k_densmat_rel(n_k_points,kdm_complex_rel,density_matrix)
       else

          if (real_eigenvectors) then
             call aims_allocate( kdm, n_basis, n_basis,                 "kdm" )
             ! dummy for the real case
             call aims_allocate( kdm_complex, 1, 1,             "kdm_complex" )
          else
             call aims_allocate( kdm_complex, n_basis, n_basis, "kdm_complex" )
             ! dummy for the complex case
             call aims_allocate( kdm, 1, 1,                             "kdm" )
          end if
          if (use_batch_permutation > 0) then
             call aims_allocate( density_matrix_temp, n_basis*(n_basis+1)/2, &
                  "density_matrix_temp" )
             density_matrix_temp = 0.d0
          end if

          ! JW: If there is need for optimization, we could just calculate
          ! those density matrix elements which are actually needed.
          ! This would scale O(N^2) instead O(N^3).

          t_read = 0.d0
          t_write = 0.d0
          i_k = 0
          do i_k_point = 1,n_k_points
             if(myid == mod(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k = i_k+1

                if(elsi_read_dm .and. do_elsi_rw) then
                   call get_timestamps(tmp,t_read1)

                   if(real_eigenvectors) then
                      call elsi_restart_lapack(i_spin,i_k_point,kdm)
                   else
                      call elsi_restart_lapack(i_spin,i_k_point,kdm_complex)
                   endif

                   call get_times(tmp,t_read1,tot_time_matrix_io,&
                        tot_clock_time_matrix_io,.true.)
                   t_read = t_read+t_read1
                else
                   call evaluate_k_densmat(kdm,kdm_complex,occ_numbers,&
                        KS_eigenvector,KS_eigenvector_complex,i_spin,i_k_point,&
                        i_k)
                endif

                if(elsi_write_dm_this .and. do_elsi_rw) then
                   call get_timestamps(tmp,t_write1)

                   write(mat_file,"(A,I2.2,A,I6.6,A)") &
                      "D_spin_",i_spin,"_kpt_",i_k_point,".csc"

                   if(real_eigenvectors) then
                      call aims_elsi_write_mat(rwh_w,trim(mat_file),kdm)
                   else
                      call aims_elsi_write_mat(rwh_w,trim(mat_file),kdm_complex)
                   endif

                   call get_times(tmp,t_write1,tot_time_matrix_io,&
                        tot_clock_time_matrix_io,.true.)
                   t_write = t_write+t_write1
                endif

                if (use_local_index) then
                   if (use_batch_permutation > 0) then
                      ! TODO: Port this to multiple k-points (will only work for
                      !       gamma-point-only calculations!)
                      call accumulate_k_densmat(density_matrix_temp,density_matrix,&
                           force_packed,kdm,kdm_complex,i_k_point)
                   else
                      call aims_stop("Sparse local matrices not supported when &
                                      &using LAPACK eigensolver.", func)
                   end if
                else
                   call accumulate_k_densmat(density_matrix_sparse,density_matrix,&
                        force_packed,kdm,kdm_complex,i_k_point)
                end if
             endif ! myid == mod(k...
          enddo ! i_k_point

       endif ! end of rel/nonrel judgement

       if (use_local_index) then
          if (use_batch_permutation > 0) then
             ! TODO: Port this to multiple k-points (will only work for
             !       gamma-point-only calculations!)
             call get_full_local_matrix(density_matrix_sparse, i_spin, &
                  density_matrix_temp)
          else
             call aims_stop("Sparse local matrices not supported when using &
                            &LAPACK eigensolver.", func)
          end if
       else
          if(packed_matrix_format /= PM_none .or. force_packed) then
             call sync_density_matrix_sparse(density_matrix_sparse)
          else
             if(flag_rel.ne.REL_x2c .and. flag_rel.ne.REL_4c_dks)then 
              ! For X2C, synchronization has been done already
             call sync_density_matrix(density_matrix)
             endif
          end if
       end if

       ! Allocatable arrays that are tracked
       if (allocated(kdm))         call aims_deallocate( kdm,                 "kdm" )
       if (allocated(kdm_complex)) call aims_deallocate( kdm_complex, "kdm_complex" )
       if (allocated(kdm_complex_rel)) call aims_deallocate( kdm_complex_rel, "kdm_complex_rel" )
       if (allocated(density_matrix_temp)) call aims_deallocate( density_matrix_temp, "density_matrix_temp" )
    end if ! use_scalapack

    ! Output restart timings
    if(elsi_read_dm .and. do_elsi_rw) then
       if(.not. elsi_read_dm_done) then
          write(info_str,"(2X,A)") "Finished reading density matrices from file"
          call localorb_info(info_str,use_unit)
          write(info_str,"(2X,A,F10.3,A)") "| Time : ",t_read," s"
          call localorb_info(info_str,use_unit)
          write(info_str,"(A)") ""
          call localorb_info(info_str,use_unit)
       else
          write(info_str,"(2X,A)") "Density matrices have been read in earlier"
          call localorb_info(info_str,use_unit)
          write(info_str,"(A)") ""
          call localorb_info(info_str,use_unit)
       endif
    endif

    if(elsi_write_dm_this .and. do_elsi_rw) then
       write(info_str,"(2X,A)") "Finished writing density matrices to file"
       call localorb_info(info_str,use_unit)
       write(info_str,"(2X,A,F10.3,A)") "| Time : ",t_write," s"
       call localorb_info(info_str,use_unit)
       write(info_str,"(A)") ""
       call localorb_info(info_str,use_unit)
    endif

  end subroutine evaluate_densmat

  !******
  !------------------------------------------------------------------------------
  !****s* density_matrix_evaluation/evaluate_k_densmat
  !  NAME
  !    evaluate_k_densmat
  !  SYNOPSIS

  subroutine evaluate_k_densmat(kdensmat, kdensmat_complex, occ_numbers, &
  &                             KS_eigenvector, KS_eigenvector_complex, &
  &                             i_spin, i_k_point, i_k)

    !  PURPOSE
    !    Calculate the k-dependend density matrix
    !    for one k and one spin only for real_eigenvectors.
    !  USES
    use aims_memory_tracking, only : aims_allocate, aims_deallocate
    use dimensions, only: n_basis, n_states, n_spin, n_k_points_task, &
        n_k_points
    use runtime_choices, only: real_eigenvectors, occupation_type
    implicit none

    !  ARGUMENTS

    real*8, intent(OUT) :: kdensmat(n_basis, n_basis)
    complex*16, intent(OUT) :: kdensmat_complex(n_basis, n_basis)
    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
    ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
    !        properly k-weighted (and are thus not the "true" occupation numbers.)
    !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
    real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    integer, intent(IN) :: i_spin, i_k_point, i_k

    !  INPUTS
    !    o occ_numbers -- occupation of states, with k-weighting already applied
    !    o KS_eigenvector, KS_eigenvector_complex -- eigencoefficients
    !    o i_spin -- spin component
    !    o i_k_point -- k-point index
    !    o i_k -- node-local k-point index
    !  OUTPUTS
    !    o kdensmat, kdensmat_complex -- density matrix of this k-point
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_state, n_max_occupied, info, i_bas
    real*8 :: occu
    real*8, allocatable :: KS_scaled(:,:)
    complex*16, allocatable :: KS_scaled_cmplx(:,:)
    logical :: use_collective
    character(*), parameter :: func = "evaluate_k_densmat"

    if (real_eigenvectors) then
       kdensmat = 0.d0
    else
       kdensmat_complex = (0.d0, 0.d0)
    end if

    if (occupation_type == 2) then
       ! dsyrk/zherk doesn't work for MP because of negative occupation numbers
       use_collective = .false.
    else
       if (count(occ_numbers(:, i_spin, i_k_point) > 0.d0) > n_states/2) then
          use_collective = .true.
       else
          use_collective = .false.
       end if
    end if

    if (use_collective) then
       ! Most states are occupied; use collective operations.
       if (real_eigenvectors) then
          call aims_allocate(KS_scaled, n_basis, n_states, "KS_scaled")
       else
          call aims_allocate(KS_scaled_cmplx, n_basis, n_states, "KS_scaled_cmplx")
       end if

       n_max_occupied = 1

       do i_state = 1, n_states, 1
          occu = occ_numbers(i_state, i_spin, i_k_point)

          if (occu > 0.d0) then
             n_max_occupied = i_state
          end if

          if (real_eigenvectors) then
            KS_scaled(:, i_state) = KS_eigenvector(:, i_state, i_spin, i_k) * sqrt(occu)
          else
            KS_scaled_cmplx(:, i_state) = KS_eigenvector_complex(:, i_state, i_spin, i_k) * sqrt(occu)
          end if
       end do

       if (real_eigenvectors) then
          call dsyrk("U", "N", n_basis, n_max_occupied, 1.d0, KS_scaled, n_basis, &
               0.d0, kdensmat, n_basis)
       else
          call zherk("U", "N", n_basis, n_max_occupied, 1.d0, KS_scaled_cmplx, n_basis, &
               0.d0, kdensmat_complex, n_basis)
       end if
    else
       ! Most states are unoccupied; use selective operations.
       do i_state = 1, n_states
          occu = occ_numbers(i_state, i_spin, i_k_point)

          if (occu > 0.d0 .or. occupation_type == 2) then
             if (real_eigenvectors) then
                call dsyr("U", n_basis, occu, &
                     KS_eigenvector(:, i_state, i_spin, i_k), 1, kdensmat, n_basis)
             else
                call zher("U", n_basis, occu, &
                     KS_eigenvector_complex(:, i_state, i_spin, i_k), 1, &
                     kdensmat_complex, n_basis)
             end if
          end if
       end do
    end if

    if (real_eigenvectors) then
       kdensmat = kdensmat + transpose(kdensmat)

       do i_bas = 1, n_basis
          kdensmat(i_bas, i_bas) = kdensmat(i_bas, i_bas)*0.5d0
       end do
    else
       kdensmat_complex = kdensmat_complex + transpose(dconjg(kdensmat_complex))

       do i_bas = 1, n_basis
          kdensmat_complex(i_bas, i_bas) = kdensmat_complex(i_bas, i_bas)*0.5d0
       end do
    end if

   if(allocated(KS_scaled)) call aims_deallocate(KS_scaled, "KS_scaled")
   if(allocated(KS_scaled_cmplx)) call aims_deallocate(KS_scaled_cmplx, "KS_scaled_cmplx")

  end subroutine evaluate_k_densmat
  !******
  !------------------------------------------------------------------------------
  !****s* density_matrix_evaluation/accumulate_k_densmat
  !  NAME
  !    accumulate_k_densmat
  !  SYNOPSIS

  subroutine accumulate_k_densmat(density_matrix_sparse, density_matrix, force_packed, &
  &                               kdm, kdm_complex, i_k_point)

    !  PURPOSE
    !    Accumulate the k-dependent density matrices to the real-space representation.
    !  USES
    use dimensions, only: n_basis, n_centers_basis_t
    use load_balancing, only: use_batch_permutation
    use mpi_tasks, only: aims_stop
    use pbc_lists, only: k_phase, index_hamiltonian, column_index_hamiltonian, &
        Cbasis_to_basis, Cbasis_to_center,  n_cells_in_hamiltonian, &
        center_to_cell
    use runtime_choices, only: real_eigenvectors, packed_matrix_format, &
        PM_index, PM_none, use_local_index
    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: density_matrix_sparse(*)
    real*8, intent(INOUT) :: density_matrix(:,:)
    logical, intent(IN) :: force_packed
    real*8, intent(IN) :: kdm(:,:)
    complex*16, intent(IN) :: kdm_complex(:,:)
    integer, intent(IN) :: i_k_point

    !  INPUTS
    !    o density_matrix{_sparse} -- real-space density matrix
    !    o force_packed -- use (lapack-)packed density_matrix even for PM_none
    !    o kdm{_cmplx} -- k-dependent density matrix for i_k_point
    !    o i_k_point -- corresponding k_point
    !  OUTPUTS
    !    o density_matrix{_sparse} -- real-space density matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    ! Add these matrix elements with their phases to the corresponding
    ! real space matrix.

    integer, allocatable :: index_b(:,:)
    real*8 :: add
    integer :: i_bas, i_bas1, i_bas2, i_basT1, i_basT2
    integer :: i_cell1, i_cell2
    integer :: i_cell
    integer :: i_index
    complex*16 :: conj_phase
    character(*), parameter :: func = "accumulate_k_densmat"

    if (use_local_index) then
       if (use_batch_permutation > 0) then
          ! TODO: Port this to multiple k-points!
          i_index = 0
          do i_bas2 = 1,  n_basis
             do i_bas1 = 1, i_bas2
               i_index = i_index + 1
               if (real_eigenvectors) then
                   add = kdm(i_bas1, i_bas2)
               else
                   add = dble(kdm_complex(i_bas1, i_bas2))
                end if
                density_matrix_sparse(i_index) = density_matrix_sparse(i_index) + add
             end do
          end do
          return
       else
          call aims_stop("Sparse local matrices not supported when using LAPACK &
                         &eigensolver.", func)
       end if
    end if

    select case(packed_matrix_format)
    case(PM_index)
       do i_bas2 = 1, n_basis
          do i_cell = 1,n_cells_in_hamiltonian-1
             conj_phase = conjg(k_phase(i_cell,i_k_point))
             if (index_hamiltonian(1,i_cell, i_bas2) > 0) then
                do i_index = index_hamiltonian(1, i_cell, i_bas2), &
                &            index_hamiltonian(2, i_cell, i_bas2)
                   i_bas1 =  column_index_hamiltonian(i_index)

                   if (real_eigenvectors) then
                      density_matrix_sparse(i_index) &
                      & = density_matrix_sparse(i_index) &
                      & + kdm(i_bas1, i_bas2) * dble(conj_phase)
                   else
                      density_matrix_sparse(i_index) &
                      & = density_matrix_sparse(i_index) &
                      & + dble(kdm_complex(i_bas1, i_bas2) * conj_phase)
                   end if

                end do
             end if
          end do
       end do

    case(PM_none)

       i_index = 0
       do i_basT2 = 1,  n_centers_basis_T
          i_bas2 = Cbasis_to_basis(i_basT2)
          i_cell2 = center_to_cell(Cbasis_to_center(i_basT2))
          do i_basT1 = 1, n_centers_basis_T
             if (force_packed .and. i_basT1 > i_basT2) cycle
             i_bas1 = Cbasis_to_basis(i_basT1)
             i_cell1 = center_to_cell(Cbasis_to_center(i_basT1))
             i_index = i_index + 1

             if (real_eigenvectors) then
                add = kdm(i_bas1, i_bas2) &
                &     * dble(k_phase(i_cell1, i_k_point)) * dble(k_phase(i_cell2, i_k_point))
             else
                add = dble(kdm_complex(i_bas1, i_bas2) &
                &     * k_phase(i_cell1, i_k_point) * dconjg(k_phase(i_cell2, i_k_point)))
             end if
             if (force_packed) then
                density_matrix_sparse(i_index) = density_matrix_sparse(i_index) + add
             else
                density_matrix(i_basT1, i_basT2) = density_matrix(i_basT1, i_basT2) + add
             end if
          end do
       end do

    case default
       call aims_stop("Invalid packing type", func)

    end select

  end subroutine accumulate_k_densmat
  !******
  !****s* density_matrix_evaluation/output_densmat
  !  NAME
  !    output_densmat
  !  SYNOPSIS
  subroutine output_densmat(KS_eigenvector,KS_eigenvector_complex,occ_numbers)
    ! PURPOSE
    !   Output the density matrix after a converged s.c.f. cycle.
    ! USES
    use aims_memory_tracking, only: aims_allocate,aims_deallocate
    use dimensions, only: n_basis, n_states, n_spin, n_k_points_task, &
        n_k_points
    use elsi_wrapper, only: rwh_w, aims_elsi_write_mat
    use localorb_io, only: use_unit, localorb_info
    use mpi_tasks, only: myid, n_tasks
    use runtime_choices, only: use_scalapack, real_eigenvectors
    use scalapack_wrapper, only: ham, ham_complex, my_k_point, &
        construct_dm_scalapack, set_full_matrix_real, set_full_matrix_complex
    use timing_core, only: get_times, get_timestamps
    use timing, only: tot_time_matrix_io, tot_clock_time_matrix_io
    implicit none

    ! ARGUMENTS
    real*8,     dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8,     dimension(n_states,n_spin,n_k_points), intent(in)  :: occ_numbers

    ! Local variables
    character(*), parameter :: func = "output_densmat"

    real*8,     dimension(:,:), allocatable :: kdm
    complex*16, dimension(:,:), allocatable :: kdm_complex

    character*100 :: mat_file
    character*150 :: info_str

    integer :: i_spin
    integer :: i_k_point
    integer :: i_k

    real*8 :: t_write
    real*8 :: t_write1
    real*8 :: tmp

    do i_spin = 1,n_spin
       if(use_scalapack) then
          call construct_dm_scalapack(occ_numbers,i_spin)

          call get_timestamps(t_write,tmp)

          write(mat_file,"(A,I2.2,A,I6.6,A)") &
             "D_spin_",i_spin,"_kpt_",my_k_point,".csc"

          if(real_eigenvectors) then
             call set_full_matrix_real(ham(:,:,i_spin))
             call aims_elsi_write_mat(rwh_w,trim(mat_file),ham(:,:,i_spin))
          else
             call set_full_matrix_complex(ham_complex(:,:,i_spin))
             call aims_elsi_write_mat(rwh_w,trim(mat_file),&
                  ham_complex(:,:,i_spin))
          endif

          call get_times(t_write,tmp,tot_time_matrix_io,&
               tot_clock_time_matrix_io)
       else ! use_scalapack
          if(real_eigenvectors) then
             call aims_allocate(kdm,n_basis,n_basis,"kdm")
             call aims_allocate(kdm_complex,1,1,"kdm_complex")
          else
             call aims_allocate(kdm_complex,n_basis,n_basis,"kdm_complex")
             call aims_allocate(kdm,1,1,"kdm")
          endif

          t_write = 0.d0
          i_k = 0

          do i_k_point = 1,n_k_points
             if(myid == mod(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k = i_k+1

                call evaluate_k_densmat(kdm,kdm_complex,occ_numbers,&
                     KS_eigenvector,KS_eigenvector_complex,i_spin,i_k_point,i_k)

                call get_timestamps(t_write1,tmp)

                write(mat_file,"(A,I2.2,A,I6.6,A)") &
                   "D_spin_",i_spin,"_kpt_",i_k_point,".csc"

                if(real_eigenvectors) then
                   call aims_elsi_write_mat(rwh_w,trim(mat_file),kdm)
                else
                   call aims_elsi_write_mat(rwh_w,trim(mat_file),kdm_complex)
                endif

                call get_times(t_write1,tmp,tot_time_matrix_io,&
                     tot_clock_time_matrix_io,.true.)
                t_write = t_write+t_write1
             endif
          enddo ! i_k_point

          call aims_deallocate(kdm,"kdm")
          call aims_deallocate(kdm_complex,"kdm_complex")
       endif ! use_scalapack

       write(info_str,"(2X,A)") "Finished writing density matrices to file"
       call localorb_info(info_str,use_unit)
       write(info_str,"(2X,A,F10.3,A)") "| Time : ",t_write," s"
       call localorb_info(info_str,use_unit)
       write(info_str,"(A)") ""
       call localorb_info(info_str,use_unit)
    enddo ! i_spin

  end subroutine output_densmat
  !******

end module density_matrix_evaluation
!******
