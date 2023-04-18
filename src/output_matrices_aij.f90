!****s* FHI-aims/output_real_matrices_aij
!  NAME
!    output_real_matrices_aij
!  SYNOPSIS

subroutine output_real_matrices_aij &
     ( overlap_matrix, hamiltonian )

!  PURPOSE
!  Subroutine output_matrices writes out all basis function properties, and the
!  previously calculated overlap and Hamiltonian matrices in a format: row, col, value.
!
!  Output functionality for complex and real versions of matrices.
!  Split into two subroutines due to variable types, and to adjust
!  formatting.
!
!  USES

  use dimensions
  use runtime_choices
  use basis
  use localorb_io
  use pbc_lists
  use mpi_tasks, only: aims_stop, check_allocation
  implicit none

!  ARGUMENTS
!  imported variables

  real*8 hamiltonian( n_hamiltonian_matrix_size,n_spin )
  real*8 overlap_matrix( n_hamiltonian_matrix_size )

!  INPUTS
!   o overlap_matrix -- overlap matrix
!   o hamiltonian -- Hamiltonian matrix
!
!  OUTPUT
!    none
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
!    Edit 2015 (harald.oberhofer@tum.de) can now also write Gamma point Hamiltonian 
!        of periodic system
!  SOURCE




!  local variables

      real*8 output_element(n_basis)

      character*18,dimension(n_spin) :: file_name

      character*100 :: info_str

      real*8,    dimension(:,:),allocatable :: output_matrix_w
      complex*16,dimension(:,:),allocatable :: output_matrix_w_complex
      real*8,    dimension(:,:,:),allocatable :: work_ham
      
      integer :: packed_matrix_format_local ! local PM format to allow for previously unpacked matrices
      complex*16 k_phase_save(n_cells) ! buffer for phase of k-point 1
      
!  counters

      integer :: i_basis, i_basis_1, i_basis_2, i_index, i_spin, i_fn
      integer :: i_cell, i_index_real, i_size, info

      
! allocate work and output arrays, mainly necessary for periodic matrix output
      allocate( output_matrix_w( n_basis*(n_basis+1)/2, n_spin), stat=info)
      call check_allocation( info, 'output_matrix_w               ')
      if( .not. real_eigenvectors) then
         allocate( output_matrix_w_complex( n_basis*(n_basis+1)/2, n_spin), stat=info)
         call check_allocation( info, 'output_matrix_w_complex       ')
      else
         allocate( output_matrix_w_complex( 1, 1), stat=info)
         call check_allocation( info, 'output_matrix_w_complex       ')
      endif
      if(packed_matrix_format == PM_none)then
         allocate(work_ham( n_centers_basis_I, n_centers_basis_I, n_spin), stat=info)
         call check_allocation( info, 'work_ham                      ') 
      else
         ! work_ham only needed for non-packed case, here: dummy only for call below
         allocate(work_ham(1, 1, 1),stat=info)
         call check_allocation( info, 'work_ham                      ') 
      end if
      
!  begin work

      write (info_str,'(2X,A)') &
        "Writing basis function properties and matrices."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
      call localorb_info(info_str,use_unit,'(A)')

!  write basis function properties

      open (50, file="basis-indices.out")

      write(info_str,*)
      call localorb_info(info_str,50,'(A)')

      write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
        "fn.", "  type  ", "at.", "n", &
        "l", "m"
      call localorb_info(info_str,50,'(A)')

      do i_basis = 1, n_basis, 1
        i_fn = basis_fn(i_basis)

        write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
          i_basis, basisfn_type(i_fn), &
          basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
          basis_m(i_basis)
        call localorb_info(info_str,50,'(A)')

      enddo

      close (50)

      packed_matrix_format_local=packed_matrix_format !local PM state as the construct_* routines already unpack the matrices
      
      k_phase_save=k_phase(:,1) ! locally store the actual phase
      k_phase(:,1)=(1.0d0,1.0d0) ! shift k-point one to Gamma point

!  now, write the overlap matrix

      open (50, file="overlap-matrix.out")

      if( n_cells > 1) then
         call construct_overlap(overlap_matrix,output_matrix_w, output_matrix_w_complex,1,work_ham)
         packed_matrix_format_local=PM_none ! matrix is already unpacked
         if( .not. real_eigenvectors) then !at the gamma point the matrices are real
            do i_size = 1,n_basis*(n_basis+1)/2
               output_matrix_w(i_size,1)=REAL(output_matrix_w_complex(i_size,1))
            enddo
         endif
         do i_size = 1,n_basis*(n_basis+1)/2
         enddo
      else
         do i_size = 1,n_basis*(n_basis+1)/2
            output_matrix_w(i_size,1)=overlap_matrix(i_size)
         enddo
      endif

      select case( packed_matrix_format_local)
      case(PM_index)
         i_cell = 1 ! cells taken care of earlier
         do i_basis_2 = 1, n_basis

            if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

               i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

               do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                  i_index_real = i_index_real + 1
                  i_basis_1 =  column_index_hamiltonian(i_index_real)

                  write(info_str,'(I6,I6,E21.10E3)') i_basis_1, i_basis_2, output_matrix_w(i_index_real,1)
                  call localorb_info(info_str,50,'(A)')

               end do
            end if
         end do

      case(PM_none)

         do i_basis_2 = 1, n_basis
            i_basis = i_basis_2*(i_basis_2-1)/2
            do i_basis_1 = 1, i_basis_2

               write(info_str,'(I6,I6,E21.10E3)') i_basis_1, i_basis_2, output_matrix_w(i_basis+i_basis_1,1)
               call localorb_info(info_str,50,'(A)')

            end do
         end do

      case default
     
         call localorb_info('Invalid packing type')
         call aims_stop

      end select

      close(50)

!  now, write the Hamiltonian matrix
     if( n_cells > 1) then
         call construct_hamiltonian(hamiltonian,output_matrix_w, output_matrix_w_complex,1,work_ham)
         packed_matrix_format_local=PM_none ! matrix is already unpacked
         if( .not. real_eigenvectors) then !at the gamma point the matrices are real
               do i_size = 1,n_basis*(n_basis+1)/2
                  do i_spin = 1, n_spin
                      output_matrix_w(i_size,i_spin)=REAL(output_matrix_w_complex(i_size,i_spin))
                  enddo
               enddo
         endif
      else
         do i_size = 1,n_basis*(n_basis+1)/2
            do i_spin = 1, n_spin
                output_matrix_w(i_size,i_spin)=hamiltonian(i_size,i_spin)
            enddo
         enddo
      endif
      if (n_spin.eq.1) then
        write(file_name(1),'(A)') "hamiltonian.out"
      else if (n_spin.eq.2) then
        write(file_name(1),'(A)') "hamiltonian_up.out"
        write(file_name(2),'(A)') "hamiltonian_dn.out"
      end if

      do i_spin = 1, n_spin, 1

         open (50, file=file_name(i_spin))

         select case( packed_matrix_format_local)
         case(PM_index)

            i_cell = 1 ! cells taken care of earlier
            do i_basis_2 = 1, n_basis

               if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                  
                  i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

                  do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                     i_index_real = i_index_real + 1
                     i_basis_1 =  column_index_hamiltonian(i_index_real)

                     write(info_str,'(I6,I6,E21.10E3)') i_basis_1, i_basis_2, &
                          output_matrix_w(i_index_real, i_spin)
                     call localorb_info(info_str,50,'(A)')

                  end do
               end if
            end do

         case(PM_none)


            do i_basis_2 = 1, n_basis
               i_basis = i_basis_2*(i_basis_2-1)/2
               do i_basis_1 = 1, i_basis_2
                  
                  write(info_str,'(I6,I6,E21.10E3)') i_basis_1, i_basis_2, &
                       output_matrix_w(i_basis+i_basis_1,i_spin)
                  call localorb_info(info_str,50,'(A)')
               
               end do
            end do

         case default
     
            call localorb_info('Invalid packing type')
            call aims_stop

         end select
         
         close(50)
         
      enddo ! spin
      
      k_phase(:,1)=k_phase_save ! restore the actual phase
      deallocate(output_matrix_w)
      deallocate(output_matrix_w_complex)
      deallocate(work_ham)


    end subroutine output_real_matrices_aij
!----------------------------------------------------------------------
!******
