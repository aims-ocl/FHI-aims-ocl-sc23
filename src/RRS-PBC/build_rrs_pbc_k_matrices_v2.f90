!****s* FHI-aims/RRS-PBC/build_rrs_pbc_k_matrices_v2()
!  NAME
!  build_rrs_pbc_k_matrices_v2()
!  SYNOPSIS

    subroutine build_rrs_pbc_k_matrices_v2(k_vector,iop)

!  PURPOSE
!  This routine generates the hamiltonian and overlap matrices associated with
!  k_vector
!
!  USES

      use localorb_io
      use runtime_choices
      use dimensions
      use physics
      use pbc_lists
      use geometry
      use numerical_utilities
      use basis
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
      real*8, dimension(3)    :: k_vector
      real*8, dimension(3)    :: k_vector_cart
      complex*16,dimension(2) :: c_para_bak
      complex*16              :: c_para
      complex*16, parameter   :: c_i = cmplx(0.0d0,1.0d0)
      integer, optional       :: iop
      integer                 :: dirt

      ! local variables
      character*132 :: rrs_pbc_info
      character*132 :: func = 'build_rrs_pbc_k_matrices()'
      integer       :: i, j, o, p, i_spin
      integer       :: i_start, j_start, o_start, p_start
      integer       :: i_n, j_n, o_n, p_n
      integer       :: i_shift, j_shift, o_shift, p_shift
      integer       :: i_abs, j_abs, o_abs, p_abs
      integer       :: i_k, j_k, j_k_start
      integer       :: c_index, k_index
      integer       :: info
      real*8,dimension(3) :: n_vec_1, n_vec_2, n_vec_3

      if (present(iop)) then
          dirt = iop
      else
          dirt = 0
      endif

      ! generate k_vector in cartesian coordinate
      k_vector_cart = k_vector(1) * rrs_pbc_recip_lattice_vector(:,1) + &
                      k_vector(2) * rrs_pbc_recip_lattice_vector(:,2) + &
                      k_vector(3) * rrs_pbc_recip_lattice_vector(:,3)

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
                      if (j_abs < i_abs) then
                          c_index = j_abs + i_abs * (i_abs-1)/2
                      else
                          c_index = i_abs + j_abs * (j_abs-1)/2
                      endif
                      k_index = j_k + i_k * (i_k-1)/2
                      ! loading gamma value for k_index
                      do i_spin = 1, rrs_pbc_n_spin, 1
                          rrs_pbc_hamiltonian_k(k_index,i_spin) = &
                              hamiltonian(c_index,i_spin)
                      enddo
                      rrs_pbc_overlap_k(k_index) = overlap_matrix(c_index)
                      !if (i .eq. 1) then
                      !    write(use_unit,*) '==========================='
                      !    write(use_unit,*) i_abs,j_abs,c_index,i_k,j_k, k_index
                      !    write(use_unit,*) real(rrs_pbc_hamiltonian_k(k_index,1)),hamiltonian(c_index,1)
                      !endif
                      ! loading the contribution of the other equal points
                      ! for k_index
                      !-- equal points associated with the j_abs basis
                      do p = 1, rrs_pbc_n_equal, 1
                          if (rrs_pbc_equal_atom(2,p,j).eq.0) cycle
                          call get_rrs_pbc_cell_vector(j,p,n_vec_1,1)
                          call get_rrs_pbc_cell_index(j,p,n_vec_2,1)
                          p_start = rrs_pbc_equal_atom(2,p,j)
                          p_shift = j_shift
                          p_abs   = p_start + p_shift - 1
                          ! to make sure the first index should be smaller than
                          ! the second index
                          if (p_abs < i_abs) then
                              c_index = p_abs + i_abs * (i_abs-1)/2
                          else
                              c_index = i_abs + p_abs * (p_abs-1)/2
                          endif
                          c_para_bak(1)  = &
                              exp(c_i * dot_product(k_vector_cart,n_vec_1))
                          c_para_bak(2)  = &
                              exp(c_i * 2.0d0*pi * dot_product(k_vector,n_vec_2))
                          c_para = c_para_bak(2)
                          if (dirt .gt. 0) then
                            if (myid .eq. 0 .and. j .eq. 1 .and. j_shift .eq. 1) then
                               !write(rrs_pbc_info,'(A,3(F8.4))') "igor debug",&
                               !    n_vec_2(1),n_vec_2(2),n_vec_2(3)
                               !call localorb_info(rrs_pbc_info,use_unit,'(A)')
                               !write(rrs_pbc_info,'(A,2(F16.8))') 'c_para_1',real(c_para_bak(1)),aimag(c_para_bak(1))
                               !call localorb_info(rrs_pbc_info)
                               !write(rrs_pbc_info,'(A,2(F16.8))') 'c_para_2',real(c_para_bak(2)),aimag(c_para_bak(2))
                               !call localorb_info(rrs_pbc_info)
                               write(use_unit,'(A,3(F8.4))') 'n_vec_2 ',n_vec_2(1),n_vec_2(2),n_vec_2(3)
                               write(use_unit,'(A,2(F16.8))') 'c_para_1',real(c_para_bak(1)),aimag(c_para_bak(1))
                               write(use_unit,'(A,2(F16.8))') 'c_para_2',real(c_para_bak(2)),aimag(c_para_bak(2))
                            endif
                          endif
                          do i_spin = 1, rrs_pbc_n_spin, 1
                              rrs_pbc_hamiltonian_k(k_index,i_spin) = &
                                  rrs_pbc_hamiltonian_k(k_index,i_spin) + &
                                  c_para * hamiltonian(c_index,i_spin)
                          enddo
                          rrs_pbc_overlap_k(k_index) = &
                              rrs_pbc_overlap_k(k_index) + &
                              c_para * overlap_matrix(c_index)
                          !if (i .eq. 1) then
                          !    write(use_unit,*) p_abs, c_index
                          !    write(use_unit,*) real(rrs_pbc_hamiltonian_k(k_index,1)),hamiltonian(c_index,1)
                          !endif
                      enddo
                  enddo
              enddo
          enddo
      enddo

      end subroutine

