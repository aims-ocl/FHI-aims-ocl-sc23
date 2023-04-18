! The following subroutines are debug subroutines that were taken out of soc_utilities,
! as I am not sure whether they still work.

!  public ::  write_soc_matrix
!  public ::  write_soc_hamil_work
!  public ::  write_unperturb_eigenval
!  public ::  write_unperturb_eigenvec
!  public ::  write_SOC_perturbed_eigenvectors
!  public ::  write_SOC_perturbed_eigenvectors_original

  subroutine write_soc_matrix(filename, soc_matrix)
    use runtime_choices
    use dimensions
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use dimensions_soc, only : n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll
    implicit none
  
    real*8, dimension(n_hamiltonian_matrix_size, 3), intent(in) :: soc_matrix
    character(len=*), intent(in) :: filename
  
    integer :: i_cell, i_basis_2, i_size
    integer :: i_basis_1, i_fn_1, i_fn_2, i_index_real
    character l_char_1, l_char_2
    character l_to_str
    character*300                  :: info_str                                    
    character(*), parameter :: func = 'write_soc_matrix'
 
    if (n_basis_soc .ne. 2*n_basis) then                                         
      write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
      call aims_stop(info_str, func)                                                 
    end if                                                                        
    if (n_basis_soc_coll .ne. 2*n_basis) then                                         
      write(info_str,'(1X,A,A)') '* Incorrect number of collinear basis elements in SOC specified.  Exiting.'
      call aims_stop(info_str, func)                                                 
    end if                                                                       
    if (mod(n_basis_soc_coll,2) .eq. 1) then                                       
      call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
    end if                                                                         
    if (n_basis_soc_ncoll .gt. 0) then                                            
      call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                   & in SOC, exiting.', func)                                        
    end if                                                                         

    if (packed_matrix_format.eq.PM_index) then
      open(99, file=filename)
      do i_cell = 1,n_cells_in_hamiltonian-1
        do i_basis_2 = 1, n_basis_soc_coll/2
          if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
            i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1
            do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
              i_index_real = i_index_real + 1
              i_basis_1 =  column_index_hamiltonian(i_index_real)
              i_fn_1 = basis_fn(i_basis_1)
              l_char_1 = l_to_str(basis_l(i_basis_1))
              i_fn_2 = basis_fn(i_basis_2)
              l_char_2 = l_to_str(basis_l(i_basis_2))
              if (i_cell.eq.1.or.abs(soc_matrix(i_index_real,1)) > 1d-8.or.abs(soc_matrix(i_index_real,2)) > 1d-8 .or. &
                   abs(soc_matrix(i_index_real,3)) > 1d-8) then
                write (99, "(I5, I5, I5, A3, I5, A10, I3, A2, I3, A5, I5, A10, I3, A2, I3, A3, F12.8, F12.8, F12.8, F12.8)" ) &
                    i_cell, i_basis_1, i_basis_2, " ( ", basis_atom(i_basis_1), basisfn_type(i_fn_1), basisfn_n(i_fn_1), l_char_1, &
                    basis_m(i_basis_1), " ) ( ", basis_atom(i_basis_2), basisfn_type(i_fn_2), basisfn_n(i_fn_2), l_char_2, &
                    basis_m(i_basis_2), " ) ", soc_matrix(i_index_real, 1),  soc_matrix(i_index_real, 2), soc_matrix(i_index_real, 3)
              end if
            end do
          end if
        end do
      end do
      close(99)
    else if (packed_matrix_format.eq.PM_none) then
      i_index_real = 0
      i_cell = 0 ! Leave in for compatability purposes, set to 0 to denote non-packed
  
      open(99, file=filename)
      do i_basis_2 = 1,n_centers_basis_I,1
        do i_basis_1 = 1,i_basis_2,1 
          i_index_real = i_index_real+1
          i_fn_1 = basis_fn(i_basis_1)
          l_char_1 = l_to_str(basis_l(i_basis_1))
          i_fn_2 = basis_fn(i_basis_2)
          l_char_2 = l_to_str(basis_l(i_basis_2))
          if (i_cell.eq.0.or.abs(soc_matrix(i_index_real,1)) > 1d-8.or.abs(soc_matrix(i_index_real,2)) > 1d-8 .or. &
               abs(soc_matrix(i_index_real,3)) > 1d-8) then
            write (99, "(I5, I5, I5, A3, I5, A10, I3, A2, I3, A5, I5, A10, I3, A2, I3, A3, F12.8, F12.8, F12.8, F12.8)" ) i_cell, &
                 i_basis_1, i_basis_2, " ( ", basis_atom(i_basis_1), basisfn_type(i_fn_1), basisfn_n(i_fn_1), l_char_1, &
                 basis_m(i_basis_1), " ) ( ", basis_atom(i_basis_2), basisfn_type(i_fn_2), basisfn_n(i_fn_2), l_char_2, &
                 basis_m(i_basis_2), " ) ", soc_matrix(i_index_real, 1),  soc_matrix(i_index_real, 2), soc_matrix(i_index_real, 3)
          end if         
        end do
      end do
      close(99)
    else
      call aims_stop("Invalid packing specified in write_soc_matrix, exiting.")
    end if
  end subroutine write_soc_matrix
  
  
  subroutine write_soc_hamil_work(filename, soc_matrix_w, soc_matrix_w_complex, SOC_Hamiltonian)
    use runtime_choices
    use dimensions
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks, only: aims_stop
    use dimensions_soc, only : n_states_soc, n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll  
    implicit none
  
    character(len=*), intent(in) :: filename
    real*8, dimension(n_basis, n_basis, 3), intent(in) :: soc_matrix_w
    complex*16, dimension(n_basis, n_basis, 3), intent(in) :: soc_matrix_w_complex
    complex*16, dimension(n_states_soc, n_states_soc), intent(in) :: SOC_Hamiltonian
  
    integer :: i_cell, i_basis_1, i_basis_2, i_state, j_state ! iterators
    integer :: i_fn_1, i_fn_2, i_index
    character l_char_1, l_char_2
    character l_to_str
    integer, dimension(n_basis, n_basis) :: index_b
    character*300                  :: info_str                                   
    character(*), parameter :: func = 'write_soc_hamil_work'
   
    if (n_basis_soc .ne. 2*n_basis) then                                         
      write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
      call aims_stop(info_str, func)                                                 
    end if                                                                       
    if (n_basis_soc_coll .ne. 2*n_basis) then                                         
      write(info_str,'(1X,A,A)') '* Incorrect number of collinear basis elements in SOC specified.  Exiting.'
      call aims_stop(info_str, func)                                                 
    end if                                                                       
    if (mod(n_basis_soc_coll,2) .eq. 1) then                                       
      call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
    end if                                                                         
                                                                                 
    if (n_basis_soc_ncoll .gt. 0) then                                            
      call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                     & in SOC, exiting.', func)                                        
    end if                                                                         
 
    ! Set up the indexing for the soc_matrix_w's
    index_b = 0
    i_index = 0
    do i_basis_2 = 1,n_basis_soc_coll/2, 1
      do i_basis_1 = 1,i_basis_2,1
        i_index = i_index + 1
        index_b(i_basis_1, i_basis_2) = i_index
        index_b(i_basis_2, i_basis_1) = i_index
      end do
    end do
  
    open(99,file=filename)
    do i_basis_1 = 1, n_basis_soc_coll/2
      i_fn_1 = basis_fn(i_basis_1)
      l_char_1 = l_to_str(basis_l(i_basis_1))
      do i_basis_2 = 1, n_basis_soc_coll/2
        i_fn_2 = basis_fn(i_basis_2)
        l_char_2 = l_to_str(basis_l(i_basis_2))
        i_index = index_b(i_basis_1, i_basis_2) 
        if (real_eigenvectors) then
          write (99, "(I5, I5, A3, I5, A10, I3, A2, I3, A5, I5, A10, I3, A2, I3, A3, F12.8, F12.8, F12.8)" ) i_basis_1, i_basis_2, &
              " ( ", basis_atom(i_basis_1), basisfn_type(i_fn_1), basisfn_n(i_fn_1), l_char_1, basis_m(i_basis_1), " ) ( ", &
              basis_atom(i_basis_2), basisfn_type(i_fn_2), basisfn_n(i_fn_2), l_char_2, basis_m(i_basis_2), " ) ", &
              soc_matrix_w(i_basis_1, i_basis_2, 1),  soc_matrix_w(i_basis_1, i_basis_2, 2), soc_matrix_w(i_basis_1, i_basis_2, 3)
        else
          write (99, "(I5, I5, A3, I5, A10, I3, A2, I3, A5, I5, A10, I3, A2, I3, A3, F12.8, F12.8, F12.8, F12.8, F12.8, F12.8)" ) &
              i_basis_1, i_basis_2, " ( ", basis_atom(i_basis_1), basisfn_type(i_fn_1), basisfn_n(i_fn_1), l_char_1, &
              basis_m(i_basis_1), " ) ( ", basis_atom(i_basis_2), basisfn_type(i_fn_2), basisfn_n(i_fn_2), l_char_2, &
              basis_m(i_basis_2), " ) ", real(soc_matrix_w_complex(i_basis_1, i_basis_2, 1)), &
              aimag(soc_matrix_w_complex(i_basis_1, i_basis_2, 1)), real(soc_matrix_w_complex(i_basis_1, i_basis_2, 2)), &
              aimag(soc_matrix_w_complex(i_basis_1, i_basis_2, 2)), real(soc_matrix_w_complex(i_basis_1, i_basis_2, 3)), &
              aimag(soc_matrix_w_complex(i_basis_1, i_basis_2, 3))
        end if
      end do
    end do
   
    do i_state = 1, n_states_soc, 1
      do j_state = 1, n_states_soc, 1
        write(99, '(I5, I5, F12.8, F12.8)') i_state, j_state, real(SOC_Hamiltonian(i_state, j_state)), &
            aimag(SOC_Hamiltonian(i_state, j_state))
      end do
    end do
  
    close(99)
  end subroutine write_soc_hamil_work
  
  subroutine write_unperturb_eigenval(file_name, k_point, KS_eigenvalue)
    use dimensions
    implicit none
  
    real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: KS_eigenvalue
    integer, intent(in) :: k_point
    character(len=*), intent(in) :: file_name
    integer :: i_state, i_spin
  
    open(99,file=file_name)
    do i_spin = 1, n_spin, 1     
      do i_state = 1, n_states, 1
        write(99,*) i_state, i_spin, KS_eigenvalue(i_state, i_spin, k_point)
      end do
    end do
    close(99)
  end subroutine write_unperturb_eigenval
  
  subroutine write_unperturb_eigenvec(file_name, k_point, KS_eigenvector, KS_eigenvector_complex)
    use runtime_choices
    use dimensions
    use pbc_lists
    use geometry
    use basis
    implicit none
    ! Have not tested which of the preceeding are actually needed.
  
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: k_point
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task), intent(in) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task), intent(in) :: KS_eigenvector_complex
  
    integer :: i_basis_1, i_state, i_spin ! iterators
    integer :: i_fn_1
    character l_char_1
    character l_to_str
  
    open(99,file=file_name)
    do i_spin = 1, n_spin, 1
      do i_state = 1, n_states, 1
        do i_basis_1 = 1, n_basis
          i_fn_1 = basis_fn(i_basis_1)
          l_char_1 = l_to_str(basis_l(i_basis_1))
          if (real_eigenvectors) then
            write (99, "(I5, I5, A3, I5, A10, I3, A2, I3, A3, F12.8)") i_state + (i_spin-1) * n_states, & 
                 i_basis_1, " ( ", basis_atom(i_basis_1),  basisfn_type(i_fn_1), basisfn_n(i_fn_1), &
                 l_char_1, basis_m(i_basis_1), " ) ", KS_eigenvector(i_basis_1, i_state, i_spin, &
                 k_point)
            write (99, "(I5, I5, A3, I5, A10, I3, A2, I3, A3, F12.8)") i_state + (i_spin-1) * n_states,  &
                  i_basis_1, " ( ", basis_atom(i_basis_1),  basisfn_type(i_fn_1), basisfn_n(i_fn_1), &
                  l_char_1, basis_m(i_basis_1), " ) ", KS_eigenvector(i_basis_1, i_state, i_spin, &
                  k_point)
          else
            write (99, "(I5, I5, A3, I5, A10, I3, A2, I3, A3, F12.8, F12.8)") i_state + (i_spin-1) * n_states, &
                  i_basis_1, " ( ", basis_atom(i_basis_1), basisfn_type(i_fn_1), basisfn_n(i_fn_1), l_char_1,   &
                  basis_m(i_basis_1), " ) ", real(KS_eigenvector_complex(i_basis_1, i_state, i_spin, k_point)),  &
                  aimag(KS_eigenvector_complex(i_basis_1, i_state, i_spin, k_point))
          end if
        end do
      end do
    end do
    close(99)
  end subroutine write_unperturb_eigenvec
  
  ! This is the original version of SOC eigenvector output based on the SR version
  ! The output is human readable, but not machine readable
  subroutine write_SOC_perturbed_eigenvectors_original(file_name, KS_eigenvector_SOC_perturbed)
    use runtime_choices
    use dimensions
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks, only : aims_stop
    use dimensions_soc, only : n_saved_states_soc, n_basis_soc, n_basis_soc_coll, &
                               n_basis_soc_ncoll
    implicit none
    ! Have not tested which of the preceeding are actually needed.
  
    character(len=*), intent(in) :: file_name
    complex*16, dimension(n_basis_soc,n_saved_states_soc), intent(in) :: KS_eigenvector_SOC_perturbed
  
    integer :: i_basis, i_state, i_state_2, i_spin, basis_offset ! iterators
    integer :: i_fn
    character l_char
    character l_to_str
    character*300 :: info_str
    character*3 :: spin_str
  
    character(*), parameter :: func = 'write_SOC_perturbed_eigenvectors_original'
  
    ! Because we are using the basis set indexing arrays from the scalar-relativistic code,
    ! we must use the same scalar-relativistic basis set
    if (n_basis_soc .ne. 2*n_basis) then
      write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
      call aims_stop(info_str, func)
    end if
    if (mod(n_basis_soc_coll,2) .eq. 1) then                                       
      call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
    end if                                                                         
                                                                                 
    if (n_basis_soc_ncoll .gt. 0) then                                            
      call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                     & in SOC, exiting.', func)                                        
    end if                                                                         
 
    open (50, FILE=file_name)
    i_state = 1
  
    ! Write out eigenstates, 4 at a time
    do while (i_state.le.(n_saved_states_soc-4))
      write(50,'(30X,5(4X,I4,4X))') &
           (i_state_2, i_state_2 = i_state, i_state+4, 1)
  
      do i_basis = 1, n_basis_soc_coll/2, 1
        do i_spin = 1, 2, 1
          if (i_spin .eq. 1) then
            basis_offset = 0
            spin_str = " up"
          else
            basis_offset = n_basis_soc_coll/2
            spin_str = " dn"
          end if 
  
          i_fn = basis_fn(i_basis)
          l_char = l_to_str(basis_l(i_basis))
  
          write(50, &
               '(I5,1X,I5,1X,A8,I3,1X,A1,1X,I3,A3,3X,5(F10.6,2X))') &
               i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
               basisfn_n(i_fn), l_char, basis_m(i_basis), spin_str, &
               (abs(KS_eigenvector_soc_perturbed(basis_offset+i_basis, i_state_2)), &
               i_state_2 = i_state, i_state+4, 1)
        end do
      end do
  
      i_state = i_state+5
      write(50,*)
    end do
  
    ! Write out the remainder
    write(50,'(30X,5(4X,I4,4X))') &
         (i_state_2, i_state_2 = i_state, n_saved_states_soc, 1)
  
    do i_basis = 1, n_basis_soc_coll/2, 1
      do i_spin = 1, 2
          if (i_spin .eq. 1) then
            spin_str = " up"
          else
            spin_str = " dn"
          end if 
        i_fn = basis_fn(i_basis)
        l_char = l_to_str(basis_l(i_basis))
  
        write(50, &
             '(I5,1X,I5,1X,A8,I3,1X,A1,1X,I3,A3,3X,5(F10.6,2X))') &
             i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
             basisfn_n(i_fn), l_char, basis_m(i_basis), spin_str, &
             (abs(KS_eigenvector_soc_perturbed(basis_offset+i_basis, i_state_2)), &
             i_state_2 = i_state, n_saved_states_soc, 1)
      end do
    end do    
  
    close(50)
  
  end subroutine write_SOC_perturbed_eigenvectors_original
