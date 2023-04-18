!****s* FHI-aims/RRS-PBC/parse_rrs_pbc_rationality()
!  NAME
!  parse_rrs_pbc_rationality()
!  SYNOPSIS

    subroutine parse_rrs_pbc_rationality()

!  PURPOSE
!  This routine checkes the rationality of the RRS-PBC scheme 
!  k_vector
!
!  USES

    use localorb_io
    use dimensions
    use physics
    use geometry
    use basis
    use runtime_choices
    use mpi_tasks
    use lapack_wrapper
    implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

    ! imported variables
    real*8, dimension(3) :: k_vector
    real*8, dimension(3) :: k_vector_cart
    complex*16           :: c_para
    complex*16, parameter:: c_i = cmplx(0.0d0,1.0d0)

    ! local variables
    character*132 :: rrs_pbc_info
    character*132 :: func = 'parse_rrs_pbc_rationality'
    integer       :: i, j, o, p, i_spin
    integer       :: i_start, j_start, o_start, p_start
    integer       :: i_n, j_n, o_n, p_n
    integer       :: i_shift, j_shift, o_shift, p_shift
    integer       :: i_abs, j_abs, o_abs, p_abs
    integer       :: i_k, j_k, j_k_start
    integer       :: c_index, k_index, m_index
    integer       :: info
    real*8,dimension(3) :: n_vec_1, n_vec_2, n_vec_3
    real*8,dimension(:,:),allocatable :: trace_hamiltonian
    real*8,dimension(:),allocatable :: trace_overlap
    real*8,dimension(3) :: max_trace
    real*8,dimension(3) :: thresh
    logical             :: cut_flag = .true.

    thresh(1) = 5.0d-2  ! for alpha-spin hamiltonian
    thresh(2) = 5.0d-2  ! for beta-spin hamiltonian
    thresh(3) = 5.0d-2  ! for overlap


    allocate(trace_hamiltonian(rrs_pbc_n_hamiltonian_k,n_spin),stat=info)
    call check_allocation(info, 'trace_hamiltonian                   ')
    allocate(trace_overlap(rrs_pbc_n_hamiltonian_k),stat=info)
    call check_allocation(info, 'trace_overlap                       ')

    i_k = 0
    j_k = 0
    ! go around all the center atoms
    do i = 1, rrs_pbc_n_center_atom, 1
        i_start = rrs_pbc_center_atom(2,i)
        i_n     = rrs_pbc_center_atom(3,i)
        ! go around all basis for this atom i
        do i_shift = 1, i_n, 1
            i_abs = i_start + i_shift - 1
            i_k   = i_k + 1
            j_k   = 0
            do j = 1, i, 1
                j_start = rrs_pbc_center_atom(2,j)
                if (j .eq. i) then
                    j_n = i_shift
                else
                    j_n = rrs_pbc_center_atom(3,j)
                endif
                do j_shift = 1, j_n, 1
                    j_k   = j_k + 1
                    j_abs = j_start + j_shift - 1
                    c_index = j_abs + i_abs * (i_abs-1)/2
                    k_index = j_k + i_k * (i_k-1)/2
                    do p = 1, rrs_pbc_n_equal, 1
                        if (rrs_pbc_equal_atom(2,p,j).eq.0) cycle
                        call get_rrs_pbc_cell_vector(j,p,n_vec_1)
                        p_start = rrs_pbc_equal_atom(2,p,j)
                        p_shift = j_shift
                        p_abs   = p_start + p_shift - 1
                        ! to make sure the first index should be smaller than
                        ! the second index
                        if (p_abs < i_abs) then
                            m_index = p_abs + i_abs * (i_abs-1)/2
                        else
                            m_index = i_abs + p_abs * (p_abs-1)/2
                        endif
                        do i_spin = 1, n_spin, 1
                                trace_hamiltonian(k_index,i_spin) = & 
                                min(abs(hamiltonian(c_index,i_spin)), &
                                    abs(hamiltonian(m_index,i_spin)))
                        enddo
                        trace_overlap(k_index) = &
                        min(abs(overlap_matrix(c_index)), &
                            abs(overlap_matrix(m_index)))
                    enddo
                enddo
            enddo
        enddo
    enddo

    !write(rrs_pbc_info,'(4X,A)') '| RRS-PBC :: Real cut for hamiltonian'
    !call localorb_info(rrs_pbc_info)
    !call output_rrs_pbc_matrix_real(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis,&
    !    rrs_pbc_n_hamiltonian_k,trace_hamiltonian)

    !write(rrs_pbc_info,'(4X,A)') '| RRS-PBC :: Real cut for overlap'
    !call localorb_info(rrs_pbc_info)
    !call output_rrs_pbc_matrix_real(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis,&
    !    rrs_pbc_n_hamiltonian_k,trace_overlap)

    max_trace(1) = 0.0d0
    max_trace(2) = 0.0d0
    max_trace(3) = 0.0d0
    do i=1, rrs_pbc_n_hamiltonian_k, 1
        do i_spin=1, n_spin, 1
            max_trace(i_spin) = max(trace_hamiltonian(i,i_spin),max_trace(i_spin))
        enddo
        max_trace(3) = max(trace_overlap(i),max_trace(3))
    enddo

    write(rrs_pbc_info,'(4X,A,3(A16))')   '|         ','H-alpha','H-beta','Overlap'
    call localorb_info(rrs_pbc_info)
    write(rrs_pbc_info,'(4X,A,3(F16.8))') '| Thresh  ',(thresh(i), i=1,3)
    call localorb_info(rrs_pbc_info)
    write(rrs_pbc_info,'(4X,A,3(F16.8))') '| Real cut',(max_trace(i), i=1,3)
    call localorb_info(rrs_pbc_info)

    do i = 1, 3, 1
        if (thresh(i) < max_trace(i)) then
            cut_flag = .false.
            write(rrs_pbc_info,'(A)') '| FCM should be enlarged to get reasonable RRS-PBC results'
            call localorb_info(rrs_pbc_info)
            !call aims_stop_coll('', func)
        endif
    enddo

    deallocate(trace_hamiltonian)
    deallocate(trace_overlap)

    ! Now test whether the overlap matrix is positive-define matrix
    k_vector = 0.0d0
    call build_rrs_pbc_k_matrices(k_vector)

    allocate(trace_hamiltonian(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis),stat=info)
    call check_allocation(info, 'trace_hamiltonian                   ')
    do i = 1, rrs_pbc_n_center_basis, 1
        do j = 1, i, 1
            c_index = j + i * (i - 1) / 2
            trace_hamiltonian(i,j) = rrs_pbc_overlap_k(c_index)
            trace_hamiltonian(j,i) = trace_hamiltonian(i,j)
        enddo
    enddo
    ! Cholesky factorization ovlp_complex = U**H * U
    call ZPOTRF('U', rrs_pbc_n_center_basis, trace_hamiltonian, &
                rrs_pbc_n_center_basis, info)
    call check_info(info, 'ZPOTRF')
    ! Calculate the inverse of U and save it in ovlp_complex
    call ZTRTRI('U', 'N', rrs_pbc_n_center_basis, trace_hamiltonian, &
                rrs_pbc_n_center_basis, info)
    call check_info(info, 'ZTRTRI')

    deallocate(trace_hamiltonian)

    contains
       subroutine check_info(info, name)
          integer info
          character*(*) name
          character*150 :: info_str
    
          if(info/=0) then
             write(info_str,*) name,' failed, info= ', info
             call aims_stop(info_str, func)
          endif
       end subroutine check_info
  end subroutine
