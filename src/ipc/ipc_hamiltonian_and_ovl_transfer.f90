! Module handling the passing of the hamiltonian and overlap matrix for each k-point represented in the k-phase list through
! the shared library. Module is currently only used and adapted for ATK. Only addition to main aims source code is single line 
! in out_plot_band.f90 calling initiate_hamiltonian_and_ovl_ipc_transfer


MODULE ipc_hamiltonian_and_ovl_transfer

     USE dimensions
     USE ipc, only: ipc_write_array, ipc_start_transaction, ipc_read_array
     USE mpi_tasks, only : check_allocation
     USE pbc_lists
     USE runtime_choices
     USE physics


     IMPLICIT NONE
     
     real*8, dimension(:,:), allocatable :: atk_k_point_list ! list of k-points read in from atk 

     PUBLIC :: initiate_hamiltonian_and_ovl_ipc_transfer

     CONTAINS

     ! Checks to see if ATK requests the hamiltonian or overlap matrix. If so, set k_phase to match k_point_list sent
     ! from ATK and calculate and transfer hamiltonian and overlap matrix for each k-point
     SUBROUTINE initiate_hamiltonian_and_ovl_ipc_transfer()
         INTEGER :: i_k_point

         if(ipc_start_transaction("HAMILTONIAN_AND_OVL_QUERY_REQUEST")) then
            CALL ipc_read_array("ATK_K_POINT_LIST", atk_k_point_list)
 
            CALL set_k_phase_from_atk_k_point_list()

            do i_k_point = 1, n_k_points
               CALL calculate_and_transfer_hamiltonian_and_ovl(i_k_point)
            end do
         endif
     END SUBROUTINE initiate_hamiltonian_and_ovl_ipc_transfer

     ! Sets k_phase to match k_point_list sent from ATK
     SUBROUTINE set_k_phase_from_atk_k_point_list()
         USE constants, only: pi
         INTEGER, DIMENSION(2) :: k_point_list_shape
         INTEGER :: i_cell_n, i_k_point

         k_point_list_shape = SHAPE(atk_k_point_list)
         n_k_points = k_point_list_shape(1)

         ! reset k_phase so that it can be reallocated according to the atk_k_point_list
         deallocate(k_phase)
         allocate(k_phase(n_cells,n_k_points))

         do i_k_point = 1, n_k_points
            do i_cell_n = 1, n_cells
               k_phase( i_cell_n, i_k_point) = exp((0.d0,2.d0)*pi*sum(atk_k_point_list(i_k_point,:)*dble(cell_index(i_cell_n,:))))  
            end do
         end do

     END SUBROUTINE set_k_phase_from_atk_k_point_list

     ! Calculates hamiltonian and overlap matrix by calling construct_hamiltonian_and_ovl and then 
     ! sending it to ATK after being parsed to a dense matrix format
     SUBROUTINE calculate_and_transfer_hamiltonian_and_ovl(i_k_point)
         real*8,    dimension(:,:),allocatable :: hamiltonian_w
         real*8,    dimension(:),  allocatable :: overlap_matrix_w
         complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
         complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

         real*8, dimension(:,:,:),allocatable :: work_ham
         real*8, dimension(:,:),allocatable :: work_ovl
         INTEGER :: i_k_point

         ! Allocate memory for arrays
         allocate(hamiltonian_w_complex   (n_basis*(n_basis+1)/2,n_spin))
         allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2))

         ! dummy allocations for the real arrays.
         ! These trigger a compiler warning for -check pointers otherwise
         ! as they are included (but never touched!) in construct_hamiltonian_and_ovl below.
         allocate(hamiltonian_w   (1,1) )
         allocate(overlap_matrix_w (1) )

         if(packed_matrix_format == PM_none)then
            allocate(work_ham(n_centers_basis_I, n_centers_basis_I, n_spin))
            allocate(work_ovl(n_centers_basis_I, n_centers_basis_I))
         else
         ! dummy only, never touched
            allocate(work_ham( 1, 1, 1))
            allocate(work_ovl( 1, 1))
         end if


         call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                                            hamiltonian_w, overlap_matrix_w, &
                                            hamiltonian_w_complex, &
                                            overlap_matrix_w_complex, &
                                            i_k_point, &
                                            work_ham, &
                                            work_ovl)

         if(ipc_start_transaction("RECEIVE_HAMILTONIAN")) then
             CALL transfer_hamiltonian_to_memory(hamiltonian_w_complex)
         endif
         if(ipc_start_transaction("RECEIVE_OVERLAP_MATRIX")) then
             CALL transfer_overlap_to_memory(overlap_matrix_w_complex)
         endif

         deallocate(hamiltonian_w_complex)
         deallocate(overlap_matrix_w_complex)
         deallocate(hamiltonian_w)
         deallocate(overlap_matrix_w)
         deallocate(work_ham)
         deallocate(work_ovl)

     END SUBROUTINE calculate_and_transfer_hamiltonian_and_ovl

     ! Adapted output_complex_hamiltonian_matrix subroutine which writes hamiltonian to file while here hamiltonian is directly
     ! transferred to a python variable in ATK
     ! 
     ! hamiltonian_w_complex 			: condensed complex hamiltonian calculated from construct_hamiltonian_and_ovl
     ! complex_hamiltonian_transfer_matrix 	: dense n_spin x n_basis x n_basis*2 matrix containing hamiltonian transferred to ATK
     SUBROUTINE transfer_hamiltonian_to_memory(hamiltonian_w_complex)
         complex*16, dimension(n_basis*(n_basis+1)/2, n_spin), intent(in) :: hamiltonian_w_complex
         complex*16, dimension(n_spin, n_basis, n_basis*2) :: complex_hamiltonian_transfer_matrix
         complex*16  :: hamiltonian_matrix_entry
         integer,dimension(:,:),allocatable :: index_b
         integer :: i_row, i_col, i_spin
         integer :: i_entry, i_index

         allocate(index_b(n_basis,n_basis),stat=i_index)
         call check_allocation(i_index, 'index_b')

         index_b = 0
         i_index = 0
         do i_col = 1,n_basis, 1
            do i_row = 1,i_col,1
               i_index = i_index + 1
               index_b(i_row, i_col) = i_index
            end do
         end do

         do i_spin = 1, n_spin
            do i_row = 1, n_basis
              do i_col = 1, n_basis
                hamiltonian_matrix_entry = (0.0,0.0)
                i_entry = 0
                if (i_row < i_col) then
                  i_entry = index_b(i_row, i_col)
                  hamiltonian_matrix_entry = hamiltonian_w_complex(i_entry,i_spin)
                  complex_hamiltonian_transfer_matrix(i_spin, i_row, i_col*2-1) = REAL( hamiltonian_matrix_entry )
                  complex_hamiltonian_transfer_matrix(i_spin, i_row, i_col*2) = AIMAG( hamiltonian_matrix_entry )
                else
                  i_entry = index_b(i_col, i_row)
                  hamiltonian_matrix_entry = hamiltonian_w_complex(i_entry,i_spin)
                  complex_hamiltonian_transfer_matrix(i_spin, i_row, i_col*2-1) = REAL( hamiltonian_matrix_entry )
                  complex_hamiltonian_transfer_matrix(i_spin, i_row, i_col*2) = AIMAG( CONJG( hamiltonian_matrix_entry ) )
                end if

              end do 
            end do 
         end do          
         CALL ipc_write_array("HAMILTONIAN_MATRIX", DBLE(complex_hamiltonian_transfer_matrix))

         deallocate(index_b)
     END SUBROUTINE transfer_hamiltonian_to_memory


     ! Adapted output_complex_overlap_matrix subroutine which writes overlap matrix to file while here overlap matrix is directly
     ! transferred to a python variable in ATK
     ! 
     ! overlap_matrix_w_complex 		: condensed complex overlap matrix calculated from construct_hamiltonian_and_ovl
     ! complex_overlap_transfer_matrix 		: dense n_basis x n_basis*2 matrix containing overlap transferred to ATK
     SUBROUTINE transfer_overlap_to_memory(overlap_matrix_w_complex)
         complex*16, dimension(n_basis*(n_basis+1)/2), intent(in) :: overlap_matrix_w_complex 
         complex*16, dimension(n_spin, n_basis, n_basis*2) :: complex_overlap_transfer_matrix
         integer,dimension(:,:),allocatable :: index_b
         complex*16  :: overlap_matrix_entry
         integer :: i_row, i_col, i_spin
         integer :: i_entry, i_index


         allocate(index_b(n_basis,n_basis),stat=i_index)
         call check_allocation(i_index, 'index_b')

         index_b = 0
         i_index = 0
         do i_col = 1,n_basis, 1
            do i_row = 1,i_col,1
               i_index = i_index + 1
               index_b(i_row, i_col) = i_index
            end do
         end do

         do i_spin = 1, n_spin
            do i_col = 1,n_basis, 1
               do i_row = 1,n_basis, 1
                  overlap_matrix_entry = (0.0,0.0)
                  i_entry = 0
                  if (i_row < i_col) then
                     i_entry = index_b(i_row, i_col)
                     overlap_matrix_entry = overlap_matrix_w_complex(i_entry)
                  else
                     i_entry = index_b(i_col, i_row)
                     overlap_matrix_entry = CONJG(overlap_matrix_w_complex(i_entry))
                  end if
                  complex_overlap_transfer_matrix(i_spin, i_row, i_col*2-1) = REAL( overlap_matrix_entry )
                  complex_overlap_transfer_matrix(i_spin, i_row, i_col*2) = AIMAG( overlap_matrix_entry )
               end do 
            end do
         end do

         CALL ipc_write_array("OVERLAP_MATRIX", DBLE(complex_overlap_transfer_matrix))

         deallocate(index_b)

     END SUBROUTINE transfer_overlap_to_memory

END MODULE ipc_hamiltonian_and_ovl_transfer

