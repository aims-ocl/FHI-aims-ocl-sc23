!****s* FHI-aims/RRS-PBC/purge_rrs_pbc_hamiltonian()
!  NAME
!  purge_rrs_pbc_hamiltonian()
!  SYNOPSIS

    subroutine purge_rrs_pbc_hamiltonian()

!  PURPOSE
!  This routine purges the hamiltonian of cluster to be periodic-like.
!  USES

    use localorb_io
    use dimensions
    use physics
    use geometry
    use numerical_utilities
    use basis
    use runtime_choices
    use mpi_tasks
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
! equal_equal_int(l_in_equal,l_equal) contains different equal points
!           l_equal is the number of different equal points
!        l_in_equal is the point number in the same equal point
! equal_equal_real(1:2,l_equal) contains the info. of different equal points
!           l_equal is the number of different equal points
!                 1 is the vector length of this point
!                 2 is the weight of this point. It should be the same number
!                 of l_in_equal, but l_in_equal is an integer, but the weight is
!                 stored in real*8 type
    integer,dimension(:,:),allocatable    :: equal_equal_int
    real*8,dimension(:,:),allocatable     :: equal_equal_real
    real*8                                :: Thresh = 1.0d-5
    real*8                                :: Zero   = 0.0d0
    integer                               :: l_equal, i_equal
    integer                               :: l_in_equal, i_in_equal
    real*8                                :: tmp_length
    logical                               :: equal_flag = .false.
    integer                               :: m, n
    integer,dimension(:),allocatable      :: equal_equal_index
    real*8,dimension(2)                   :: tmp_hamiltonian
    real*8                                :: tmp_overlap

    ! local variables
    character*132         :: rrs_pbc_info
    character*132         :: func = 'purge_rrs_pbc_hamiltonian()'
    integer               :: i, j, o, p, i_spin
    integer               :: i_start, j_start, o_start, p_start
    integer               :: i_n, j_n, o_n, p_n
    integer               :: i_shift, j_shift, o_shift, p_shift
    integer               :: i_abs, j_abs, o_abs, p_abs
    integer               :: i_k, j_k, j_k_start
    integer               :: c_index, k_index
    integer               :: info
    real*8,dimension(3)   :: n_vec_1, n_vec_2, n_vec_3

    !if (allocated(equal_equal_int)) deallocate(equal_equal_int)
    !if (allocated(equal_equal_real)) deallocate(equal_equal_real)
    allocate(equal_equal_int(rrs_pbc_n_equal,rrs_pbc_n_equal))
    allocate(equal_equal_real(2,rrs_pbc_n_equal))
    write(use_unit,*) rrs_pbc_equal_atom

    ! go around all the center atoms
    do i = 1, rrs_pbc_n_center_atom, 1
        i_start = rrs_pbc_center_atom(2,i)
        i_n     = rrs_pbc_center_atom(3,i)
        do j = 1, i, 1
            j_start = rrs_pbc_center_atom(2,j)
            ! purge the weight of equal points
            do i_equal = 1, rrs_pbc_n_equal, 1
                do i_in_equal = 1, rrs_pbc_n_equal, 1
                    equal_equal_int(i_in_equal,i_equal) = 0
                    equal_equal_real(1,i_equal) = Zero
                    equal_equal_real(2,i_equal) = Zero
                enddo
            enddo
            l_equal = 0
            l_in_equal = 0
            ! go around all equal atoms
            do p = 1, rrs_pbc_n_equal, 1
                if (rrs_pbc_equal_atom(2,p,j).eq.0) cycle
                call get_rrs_pbc_cell_vector(j,p,n_vec_1)
                tmp_length = sqrt(dot_product(n_vec_1,n_vec_1))
                equal_flag = .false.
                do i_equal = 1, l_equal, 1
                    if (dabs(tmp_length - equal_equal_real(1,i_equal)) < Thresh) then
                        equal_flag = .true.
                        equal_equal_real(2,i_equal) = equal_equal_real(2,i_equal) + 1.0d0
                        equal_equal_int(int(equal_equal_real(2,i_equal)),i_equal) = p
                        exit
                    endif
                enddo
                if (.not. equal_flag) then
                    l_equal = l_equal + 1
                    equal_equal_int(1, l_equal)  = p
                    equal_equal_real(1, l_equal) = tmp_length
                    equal_equal_real(2, l_equal) = 1.0d0
                endif
            enddo
            write(rrs_pbc_info,'(4X,A,I5)') &
                '| Igor debug : Consider the equal atom of the specie ', &
                rrs_pbc_center_atom(1,j)
            call localorb_info(rrs_pbc_info)
            write(rrs_pbc_info,'(4X,A,2(I5))') &
                '| Igor debug : the first and second equal atoms which have same lattice vector length', &
                rrs_pbc_equal_atom(1,equal_equal_int(1,1),j),rrs_pbc_equal_atom(1,equal_equal_int(2,1),j)
            call localorb_info(rrs_pbc_info)
            ! go around all basis for this atom i
            do i_shift = 1, i_n, 1
                i_abs = i_start + i_shift - 1
                if (j .eq. i) then
                    j_n = i_shift
                else
                    j_n = rrs_pbc_center_atom(3,j)
                endif
                do j_shift = 1, j_n, 1
                    do i_equal = 1, l_equal, 1
                        l_in_equal = int(equal_equal_real(2,i_equal))
                        if (allocated(equal_equal_index)) deallocate(equal_equal_index)
                        allocate(equal_equal_index(l_in_equal))
                        tmp_hamiltonian = Zero
                        tmp_overlap = Zero
                        do i_in_equal = 1, l_in_equal, 1
                            p = equal_equal_int(i_in_equal, i_equal)
                            p_start = rrs_pbc_equal_atom(2,p,j)
                            p_shift = j_shift
                            p_abs   = p_start + p_shift - 1
                            ! to make sure the first index should be smaller than
                            ! the second index
                            !write(use_unit,*) "Igor Debug 0", j_n,p_start, p_shift, p_abs
                            if (p_abs < i_abs) then
                                c_index = p_abs + i_abs * (i_abs-1)/2
                            else
                                c_index = i_abs + p_abs * (p_abs-1)/2
                            endif
                            equal_equal_index(i_in_equal) = c_index
                            !write(use_unit,*) "Igor Debug 1", i_abs, c_index, equal_equal_real(2,i_equal)
                            do i_spin = 1, n_spin, 1
                                tmp_hamiltonian(i_spin) = tmp_hamiltonian(i_spin) + &
                                    dabs(hamiltonian(c_index,i_spin)/equal_equal_real(2,i_equal))
                            enddo
                            tmp_overlap = tmp_overlap + dabs(overlap_matrix(c_index)/equal_equal_real(2,i_equal))
                            write(rrs_pbc_info,'(4X,A,2(I5),2(F16.8))') &
                                '| Igor debug :', p_abs, i_abs, overlap_matrix(c_index),tmp_overlap
                            call localorb_info(rrs_pbc_info)
                        enddo
                        !do i_in_equal = 1, l_in_equal, 1
                        !    c_index = equal_equal_index(i_in_equal)
                        !    do i_spin = 1, n_spin, 1
                        !        if (hamiltonian(c_index,i_spin)*tmp_hamiltonian(i_spin)<0.0) then
                        !            hamiltonian(c_index,i_spin) = -tmp_hamiltonian(i_spin)
                        !        else
                        !            hamiltonian(c_index,i_spin) = tmp_hamiltonian(i_spin)
                        !        endif
                        !    enddo
                        !    if (overlap_matrix(c_index)*tmp_overlap<0.0) then
                        !        overlap_matrix(c_index) = -tmp_overlap
                        !    else
                        !        overlap_matrix(c_index) = tmp_overlap
                        !    endif
                        !enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    end subroutine

