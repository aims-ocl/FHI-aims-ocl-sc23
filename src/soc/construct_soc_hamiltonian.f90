!****f* FHI-aims/construct_soc_hamiltonian
!*  NAME
!*    construct_soc_hamiltonian
!*  SYNOPSIS
subroutine construct_SOC_hamiltonian(ld_soc_matrix, soc_matrix, &
                                     n_eigenvec_rows, n_eigenvec_cols, &
                                     eigenvec, eigenvec_complex, &
                                     k_point_global, &
                                     n_soc_rows, n_soc_cols, soc_ham)
!*  PURPOSE
!*    Calculates the SOC Hamiltonian matrix elements between unperturbed
!*    eigenstates to be used for the second-variational step
!*  USES
  use pbc_lists, only: n_cells_in_hamiltonian, index_hamiltonian, &
      column_index_hamiltonian, k_phase, Cbasis_to_basis, Cbasis_to_center, &
      center_to_cell
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  use dimensions, only: n_basis, n_centers_basis_I, n_spin
  use runtime_choices, only: PM_index, PM_none, packed_matrix_format, &
      real_eigenvectors, use_local_index, use_load_balancing, use_scalapack
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use scalapack_wrapper, only: sc_desc, dlen_, &
      construct_hamiltonian_like_matrix_scalapack, mxld, mxcol, l_row, l_col, &
      my_scalapack_comm_all, set_sparse_local_matrix_scalapack_generic, &
      set_full_local_matrix_scalapack_generic
  use synchronize_mpi_basic, only: sync_vector_complex
  use load_balancing, only: batch_perm, n_bp_integ, init_comm_full_local_matrix
  use scalapack_soc, only: mxld_soc, mxcol_soc, sc_desc_soc
  use dimensions_soc, only: n_states_sr, n_states_soc, sr_state_start, &
      n_core_states_omit_from_soc
  implicit none
!*  ARGUMENTS
  integer, intent(in) :: ld_soc_matrix
  real*8,  intent(in) :: soc_matrix(ld_soc_matrix, 3)
  integer, intent(in) :: n_eigenvec_rows
  integer, intent(in) :: n_eigenvec_cols
  real*8,  intent(in) :: eigenvec(n_eigenvec_rows, n_eigenvec_cols, n_spin)
  complex*16, intent(in) :: &
       eigenvec_complex(n_eigenvec_rows, n_eigenvec_cols, n_spin)
  integer, intent(in) :: k_point_global
  integer, intent(in) :: n_soc_rows
  integer, intent(in) :: n_soc_cols
  complex*16, intent(out) :: soc_ham(n_soc_rows, n_soc_cols)
!*  INPUTS
!*    o ld_soc_matrix    - Leading dimension for soc_matrix.  For all code paths except local indexing + load balancing,
!*                         n_hamiltonian_matrix_size.  For local indexing + load balancing,
!*                         batch_perm(n_bp_integ)%n_local_matrix_size
!*    o soc_matrix       - the SOC interaction strength between basis elements for a given spin coupling/Pauli matrix
!*    o n_eigenvec_rows  - Number of rows in eigenvec{_complex}.  For LAPACK, should be n_basis.
!*                         For ScaLAPACK, should be mxld.
!*    o n_eigenvec_cols  - Number of rows in eigenvec{_complex}.  For LAPACK, should be n_states.
!*                         For ScaLAPACK, should be mxcol.
!*    o eigenvec         - the unperturbed KS eigenvectors for the real eigenvectors case.
!*    o eigenvec_complex - the unperturbed KS eigenvectors for the complex eigenvectors case.
!*    o k_point_global   - the k-point to be calculated as global index (for indexing in k_phase)
!*    o n_soc_rows       - Number of rows in soc_ham.  For LAPACK, should be n_states_soc.  For ScaLAPACK, should be mxld_soc.
!*    o n_soc_cols       - Number of columns in soc_ham.  For LAPACK, should be n_states_soc.  For ScaLAPACK, should be mxcol_soc.
!*  OUTPUT
!*    o soc_ham          - The SOC Hamiltonian matrix elements between unperturbed eigenstates, in either global/LAPACK
!*                         form or local/ScaLAPACK 2D block cyclic form
!*  AUTHOR
!*    William Huhn
!*  NOTES
!*    The first part of this subroutine is a modification of construct_hamiltonian()
!*
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*
!*    This subroutine implements step 4 (non-periodic) and steps 4-5 (periodic)
!*    in Section III.3 from said paper.  Note that, for the construction of the
!*    Hamiltonian entering into the eigensolver in FHI-aims, there is no
!*    distinction in aims between periodic and non-periodic systems; a
!*    non-periodic system entering into the subroutine is treated as a periodic
!*    system with a single unit cell in the supercell and a single k-point of
!*    unity.  The distinction between the two cases in the paper as presented
!*    in the paper was done to aid readers not familiar with periodic boundary
!*    conditions, and the periodic case presented in the paper is the
!*    relevant code path for both periodic and non-periodic systems.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!* SOURCE
  character*300 :: info_str
  integer,dimension(:,:),allocatable :: index_b
  integer,dimension(:,:),allocatable :: index_b_old ! Index over Cbasis x Cbasis space
  integer,dimension(:,:),allocatable :: index_b_new ! Index over basis x basis space
  complex*16, dimension(:,:),allocatable :: work
  real*8, dimension(:,:,:), allocatable :: work_soc_matrix

  complex*16, dimension(:,:), allocatable :: temp_complex
  complex*16, dimension(:,:,:), allocatable :: temp_eigenvec_complex

  real*8 :: soc_matrix_w (1,1)
  complex*16, dimension(:,:), allocatable :: soc_matrix_w_complex ! Like soc_matrix, but between basis elements
                                                                  ! in the Bloch basis, complex case
  integer, parameter :: n_pauli = 3  ! The index corresponding to sigma_x, sigma_y, sigma_z matrices

  integer:: i_cell, i_basis, i_index_real, i_size
  integer:: i_index_new, i_index_old, i_index, i_basis_2, i_basis_1

  integer:: i_pauli
  integer :: i_state, j_state
  integer :: up_index, down_index, n_spin_save
  logical :: real_eigenvectors_save

  integer :: info

  ! Intermediate storage arrays whose size vary based on whether we are working
  ! with local or global arrays
  if (use_scalapack) then
    call aims_allocate(soc_matrix_w_complex, mxld, mxcol, "soc_matrix_w_complex")
    ! Because n_states_sr < n_states and leading dimension is the same as the
    ! SR case, we can use same communicator and matrix layout as in SR case
    ! for this matrix
    call aims_allocate(temp_complex, mxld, mxcol, "temp_complex")
  else
    call aims_allocate(soc_matrix_w_complex, n_basis, n_basis,"soc_matrix_w_complex")
    call aims_allocate(temp_complex, n_basis, n_states_sr, "temp_complex")
  end if

  ! Create a complex copy of the eigenvectors, if they are not already complex
  ! This shouldn't be necessary.  There should be a way to play with
  ! real/imaginary parts individually to work with real eigenvectors and
  ! eliminate this matrix copy.
  if (real_eigenvectors) then
    call aims_allocate(temp_eigenvec_complex, n_eigenvec_rows, n_eigenvec_cols,&
         n_spin, "temp_eigenvec_complex")
    temp_eigenvec_complex = dcmplx(eigenvec)
  else
    call aims_allocate(temp_eigenvec_complex, 1, 1, 1, "temp_cmplx_eienvec")
  end if

  soc_ham = (0.0d0, 0.0d0)
  temp_complex = (0.0d0, 0.0d0)

  if (use_scalapack.and.use_local_index) then
    if (use_load_balancing) then
      write(info_str,'(2X,A)')  'Constructing SOC Hamiltonian using local indexing + load balancing'
      call localorb_info( info_str )
    else
      write(info_str,'(2X,A)')  'Constructing SOC Hamiltonian using local indexing'
      call localorb_info( info_str )
    end if
  end if

  do i_pauli = 1, n_pauli

    soc_matrix_w_complex = cmplx(0.0d0, 0.0d0)

    ! Convert the real-space SOC matrix into the Bloch-space SOC matrix
    if (use_scalapack) then
      ! We need to trick the ScaLAPACK Hamiltonian construction into always
      ! taking the complex code path (see below for why) and bypass the
      ! spin-polarization infrastructure
      ! Because soc_matrix has an identical form to the Hamiltonian, as the spin
      ! dependence has been factored out, we can use the same ScaLAPACK
      ! descriptor AND use_local_index scheme for it.
      ! Here we do not use the SOC-to-SR environment conversion subroutines,
      ! because soc_matrix is defined on the spacial portion of the basis set,
      ! which is identical to the original scalar-relativistic basis set

      real_eigenvectors_save = real_eigenvectors
      real_eigenvectors = .false.
      if (use_local_index) then
        if (use_load_balancing) then
          call init_comm_full_local_matrix( &
               batch_perm(n_bp_integ)%n_basis_local, &
               batch_perm(n_bp_integ)%i_basis_local )
          call set_full_local_matrix_scalapack_generic(soc_matrix(1,i_pauli), &
               soc_matrix_w, soc_matrix_w_complex)
        else
          call set_sparse_local_matrix_scalapack_generic(soc_matrix(1,i_pauli),&
               soc_matrix_w, soc_matrix_w_complex )
        end if
      else
        n_spin_save = n_spin
        n_spin = 1
        call construct_hamiltonian_like_matrix_scalapack(soc_matrix(1,i_pauli),&
             soc_matrix_w, soc_matrix_w_complex )
        n_spin = n_spin_save
      end if
      real_eigenvectors = real_eigenvectors_save

      soc_matrix_w_complex = (0.0d0,-1.0d0) * soc_matrix_w_complex
    else ! The LAPACK case
      select case( packed_matrix_format )
        case(PM_index)
          ! Here, we add the k-point dependence to soc_matrix and convert to
          ! generalized Bloch coordinates.  This part is virtually identical to
          ! construct_hamiltonian.

          if (allocated(index_b)) then
            call aims_deallocate(index_b, "index_b")
          end if
          call aims_allocate(index_b, 1, 1, "index_b")
          if (allocated(work_soc_matrix)) then
            call aims_deallocate(work_soc_matrix,  "work_soc_matrix")
          end if
          call aims_allocate(work_soc_matrix, 1, 1, 1, "work_soc_matrix")

          i_index_real = 0

          ! Where we convert from the computational basis to generalized Bloch
          ! basis. i_basis_2 runs over all elements in the supercell, but
          ! i_basis_1 runs over all elements in the unit cell.  So many
          ! i_basis_2's will contribute to a single element of soc_matrix_w*
          ! with appropriate phase factor, but only a single i_basis_1 element
          ! will, as all elements of soc_matrix_w* map to a basis element in the
          ! unit cell.
          do i_cell = 1,n_cells_in_hamiltonian-1
            do i_basis_2 = 1, n_basis
              if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1
                do i_size = index_hamiltonian(1,i_cell, i_basis_2), &
                     index_hamiltonian(2,i_cell, i_basis_2)
                  i_index_real = i_index_real + 1
                  i_basis_1 =  column_index_hamiltonian(i_index_real)
                  soc_matrix_w_complex(i_basis_1,i_basis_2) = &
                    soc_matrix_w_complex(i_basis_1,i_basis_2) &
                    + (0.d0,-1.d0) * k_phase( i_cell,k_point_global)  &
                    * soc_matrix(  i_index_real , i_pauli)
                end do
              end if
            end do
          end do

        case(PM_none)
          ! Note:  index_b is identical to index_b_new from the old routine, and
          !        is identicial to index_b from the packed case, so I've
          !        renamed it here.  Should move some allocations and
          !        initialization outside the cases, but I'm feeling lazy.

          if (allocated(index_b_old)) then
            call aims_deallocate(index_b_old, &
                 "index_b_old")
          end if
          call aims_allocate(index_b_old, n_centers_basis_I, n_centers_basis_I,&
               "index_b_old")
          if (allocated(index_b)) then
            call aims_deallocate(index_b, "index_b")
          end if
          call aims_allocate(index_b, n_basis, n_basis, "index_b")

          index_b_old = 0
          index_b = 0
          i_index = 0

          ! The indexing for soc_matrix_w*
          do i_basis_2 = 1,n_basis, 1
            do i_basis_1 = 1,i_basis_2,1
              i_index = i_index + 1
              index_b(i_basis_1, i_basis_2) = i_index
              index_b(i_basis_2, i_basis_1) = i_index
            end do
          end do

          i_index = 0
          ! The indexing for soc_matrix
          do i_basis_2 = 1, n_centers_basis_I, 1
            do i_basis_1 = 1,i_basis_2,1
              i_index = i_index + 1
              index_b_old(i_basis_1, i_basis_2) = i_index
            end do
          end do

          if (allocated(work_soc_matrix)) then
            call aims_deallocate(work_soc_matrix, &
                 "work_soc_matrix")
          end if
          call aims_allocate(work_soc_matrix, 1, 1, 1, "work_soc_matrix")
        ! For real eigenvectors with no packing, either index may correspond to
        ! a basis element in the super cell, so we need to consider both
        ! i_basis_1 and i_basis_2 having a possible phase factor
          do i_basis_2 = 1,n_centers_basis_I,1
            do i_basis_1 = 1,i_basis_2,1
              i_index_old =  max(index_b_old(i_basis_2, i_basis_1), &
                   index_b_old( i_basis_1, i_basis_2))
              i_index_new = index_b(Cbasis_to_basis(i_basis_1), &
                   Cbasis_to_basis(i_basis_2))
              if(i_index_new /= 0 .and.  i_index_old /= 0) then
                ! WPH:  This functionality has not been tested extensively!
                soc_matrix_w_complex(i_basis_1,i_basis_2) =  &
                     soc_matrix_w_complex(i_basis_1,i_basis_2)   &
                     + (0.0,-1.0d0) * &
                     k_phase(center_to_cell(Cbasis_to_center(i_basis_2)), &
                     k_point_global) * &
                     k_phase(center_to_cell(Cbasis_to_center(i_basis_1)), &
                     k_point_global) * &
                     soc_matrix( i_index_old,i_pauli)
                continue
              end if
            end do
          end do
        case default
           call localorb_info('Invalid packing type')
           call aims_stop
      end select
    end if ! use_scalapack

    ! The normal Hamiltonian computation ends here with the Hamiltonian in the
    ! Bloch basis set.  But for the second-variational step we need matrix
    ! elements between unperturbed eigenstates of the system, SOC_Hamiltonian,
    ! so we compute these now.

    ! These are indices for the spin indices of the eigenvectors, used to handle
    ! both "spin none" and "spin collinear" runs in the same framework.  They
    ! have the side benefit of making the code a little easier to convert back
    ! to physical equations.
    if (n_spin.eq.1) then
      up_index = 1
      down_index = 1
    else
      up_index = 1
      down_index = 2
    end if

    ! Calculate <Psi'|V_SOC|Psi> to generate the matrix elements for the full
    ! SOC operator
    if (use_scalapack) then
      if ( (n_soc_rows.ne.mxld_soc) .or. (n_soc_cols.ne.mxcol_soc) ) then
        call aims_stop("Array indexing error in ScaLAPACK branch of construct_soc_hamiltonian!  Exiting.")
      end if
      if (real_eigenvectors) then
        if (i_pauli .eq. 1) then
          ! spin up-spin down coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,n_states_sr+1,sc_desc_soc)
          ! spin down-spin up coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,1,sc_desc_soc)
        else if (i_pauli .eq. 2) then
          ! spin up-spin down coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (0.d0,-1.d0), &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,n_states_sr+1,sc_desc_soc )
          ! spin down-spin up coupling
          call pzhemm("L", "U",  &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (0.d0,1.d0), &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,1,sc_desc_soc )
        else if (i_pauli .eq. 3) then
          ! spin up-spin up coupling
          call pzhemm('L', 'U', &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm('C', 'N', &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      temp_eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,1,sc_desc_soc )
          ! spin down-spin down coupling
          call pzhemm('L', 'U', &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm('C', 'N', &
                      n_states_sr, n_states_sr, n_basis, &
                      (-1.d0,0.d0), &
                      temp_eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,n_states_sr+1,sc_desc_soc )
        end if
      else
        if (i_pauli .eq. 1) then
          ! spin up-spin down coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,n_states_sr+1,sc_desc_soc )
          ! spin down-spin up coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,1,sc_desc_soc )
        else if (i_pauli .eq. 2) then
          ! spin up-spin down coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (0.d0,-1.d0), &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,n_states_sr+1,sc_desc_soc )
          ! spin down-spin up coupling
          call pzhemm("L", "U", &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm("C", "N", &
                      n_states_sr, n_states_sr, n_basis, &
                      (0.d0,1.d0), &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,1,sc_desc_soc )
        else if (i_pauli .eq. 3) then
          ! spin up-spin up coupling
          call pzhemm('L', 'U', &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm('C', 'N', &
                      n_states_sr, n_states_sr, n_basis, &
                      (1.d0,0.d0), &
                      eigenvec_complex(:,:,up_index), 1,sr_state_start,sc_desc,&
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, 1,1,sc_desc_soc )
          ! spin down-spin down coupling
          call pzhemm('L', 'U', &
                      n_basis, n_states_sr, &
                      (1.d0,0.d0), &
                      soc_matrix_w_complex, 1,1,sc_desc, &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      (0.d0,0.d0), &
                      temp_complex, 1,1,sc_desc)
          call pzgemm('C', 'N', &
                      n_states_sr, n_states_sr, n_basis, &
                      (-1.d0,0.d0), &
                      eigenvec_complex(:,:,down_index), 1,sr_state_start,sc_desc, &
                      temp_complex, 1,1,sc_desc, &
                      (1.d0,0.d0), &
                      soc_ham, n_states_sr+1,n_states_sr+1,sc_desc_soc )
        end if
      end if ! real_eigenvectors
    else ! LAPACK
      if ( (n_soc_rows.ne.n_states_soc) .or. (n_soc_cols.ne.n_states_soc) ) then
        call aims_stop("Array indexing error in LAPACK branch of construct_soc_hamiltonian!  Exiting.")
      end if
      if (real_eigenvectors) then
        if (i_pauli .eq. 1) then
          ! spin up-spin down coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,n_states_sr+1), n_states_soc)
          ! spin down-spin up coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,1), n_states_soc)
        else if (i_pauli .eq. 2) then
          ! spin up-spin down coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (0.d0,-1.d0), &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,n_states_sr+1), n_states_soc)
          ! spin down-spin up coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (0.d0,1.d0), &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,1), n_states_soc)
        else if (i_pauli .eq. 3) then
          ! spin up-spin up coupling
          call zhemm('L', 'U', &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm('C', 'N', &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     temp_eigenvec_complex(1,sr_state_start,up_index), n_basis,&
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,1), n_states_soc)
          ! spin down-spin down coupling
          call zhemm('L', 'U', &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm('C', 'N', &
                     n_states_sr, n_states_sr, n_basis, &
                     (-1.d0,0.d0), &
                     temp_eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,n_states_sr+1), n_states_soc)
        end if
      else
        if (i_pauli .eq. 1) then
          ! spin up-spin down coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,n_states_sr+1), n_states_soc)
          ! spin down-spin up coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,1), n_states_soc)
        else if (i_pauli .eq. 2) then
          ! spin up-spin down coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (0.d0,-1.d0), &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,n_states_sr+1), n_states_soc)
          ! spin down-spin up coupling
          call zhemm("L", "U", &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm("C", "N", &
                     n_states_sr, n_states_sr, n_basis, &
                     (0.d0,1.d0), &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,1), n_states_soc)
        else if (i_pauli .eq. 3) then
          ! spin up-spin up coupling
          call zhemm('L', 'U', &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm('C', 'N', &
                     n_states_sr, n_states_sr, n_basis, &
                     (1.d0,0.d0), &
                     eigenvec_complex(1,sr_state_start,up_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(1,1), n_states_soc)
          ! spin down-spin down coupling
          call zhemm('L', 'U', &
                     n_basis, n_states_sr, &
                     (1.d0,0.d0), &
                     soc_matrix_w_complex, n_basis, &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     (0.d0,0.d0), &
                     temp_complex, n_basis)
          call zgemm('C', 'N', &
                     n_states_sr, n_states_sr, n_basis, &
                     (-1.d0,0.d0), &
                     eigenvec_complex(1,sr_state_start,down_index), n_basis, &
                     temp_complex, n_basis, &
                     (1.d0,0.d0), &
                     soc_ham(n_states_sr+1,n_states_sr+1), n_states_soc)
        end if
      end if ! real_eigenvectors
    end if ! use_scalapack
  end do ! i_pauli

  ! Deallocate tracked arrays
  if (use_scalapack) then
    call aims_deallocate(soc_matrix_w_complex, "soc_matrix_w_complex")
    call aims_deallocate(temp_complex, "temp_complex")
  else
    call aims_deallocate(soc_matrix_w_complex, "soc_matrix_w_complex")
    call aims_deallocate(temp_complex, "temp_complex")
  end if

  if (.not.use_scalapack) then
    select case( packed_matrix_format )
      case(PM_index)
        call aims_deallocate(index_b, "index_b")
        call aims_deallocate(work_soc_matrix, "work_soc_matrix")
      case(PM_none)
        call aims_deallocate(index_b_old, "index_b_old")
        call aims_deallocate(index_b, "index_b")
        call aims_deallocate(work_soc_matrix, "work_soc_matrix")
      case default
         call localorb_info('Invalid packing type')
         call aims_stop
    end select
  end if ! use_scalapack

  ! Solely to make output prettier
  if (use_scalapack.and.use_local_index) then
    write(info_str,'(2X,A)') ''
    call localorb_info( info_str )
  end if

  call aims_deallocate(temp_eigenvec_complex, "temp_eigenvec_complex")
end subroutine construct_soc_hamiltonian
!******
