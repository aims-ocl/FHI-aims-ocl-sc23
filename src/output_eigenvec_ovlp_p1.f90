

      subroutine out_eigenvec_ovlp_p1 ()

!  PURPOSE
!  Output eigenvector and overlap matrix for periodic systems
!  By Yong Xu, Dec. 07, 2011

      use dimensions
      use constants
      use mpi_tasks
      use localorb_io
      use basis
      use geometry
      use species_data
      use runtime_choices
      use pbc_lists
      use synchronize_mpi
      use density_matrix_evaluation
      use physics
      use scalapack_wrapper
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
!    Release version, FHI-aims (2008).
!  SOURCE


      ! Local variables

      character*40 :: filename_eigenvec, filename_ovlp


      ! begin work


      write(filename_eigenvec,'(A12)') 'eigenvec.out'

      call output_eigenvalue_eigenvector_p1  & 
           ( KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, chemical_potential, filename_eigenvec )
      

      write(filename_ovlp,'(A12)') 'ovlp_mat.out'

      call output_overlap_matrix_p1 ( overlap_matrix, filename_ovlp )


    end subroutine out_eigenvec_ovlp_p1




subroutine output_eigenvalue_eigenvector_p1 &
         ( KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, chemical_potential, filename )

!  PURPOSE
!  The subroutine outputs eigenvalues&eigenvectors for periodic systems.
!  If using scalapack, the code will open the option of collecting eigenvectors.
!  See below for the detailed information about the output file.
!******************************************************************************************************

!Output eigenvectors of one k_point into one unformatted file, which has format:
!line 1:  n_basis, n_states, n_spin, n_k_points, n_atoms
!line 2:  basis_atom(1:n_basis)
!line 3 to (n_states * n_spin + 2):  
!do i_spin = 1, n_spin, 1
!   do i_state = 1, n_states, 1
!         write: eigenvalue-chemical_potential, eigenvector(1:n_basis, i_states, i_spin, i_k_point)
!   end do
!end do

!******************************************************************************************************
!
!  USES

  use dimensions
  use constants
  use mpi_tasks
  use localorb_io
  use basis
  use geometry
  use species_data
  use runtime_choices
  use pbc_lists
  use synchronize_mpi
  use density_matrix_evaluation
  use scalapack_wrapper
  implicit none


!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points)               :: KS_eigenvalue 
  real*8, dimension(n_basis, n_states, n_spin, n_k_points_task) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task) :: KS_eigenvector_complex
  real*8:: chemical_potential
  character*40 :: filename, filename_new


!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o KS_eigenvector -- Kohn-Sham eigenvectors if real format
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors if complex format
!   o chemical_potential -- chemical potential
!   o filename -- file where results are printed
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
!  SOURCE
!

  ! local variables
  character*100 :: info_str
  character*3 :: c_num

  ! counters
  integer :: i_spin, i_state, i
  integer :: i_k, i_k_point
  integer :: fileID, IRECL, iline
  integer, allocatable ::  myid_tmp(:)
  integer :: LEN, itmp
  
  INQUIRE (IOLENGTH=LEN) itmp


  ! begin work

  write(info_str,'(2X,A)') &
       ' '
  call localorb_info(info_str)
  write(info_str,'(2X,A,A)') &
       'Output eigenvectors into unformatted file: ', filename
  call localorb_info(info_str)

  if(n_periodic == 0 .and. packed_matrix_format==PM_none)then
     write(use_unit,*) 'Error: output_eigenvalue_eigenvector_p1 supports only periodic systems.'
     return
  
  else

     if(packed_matrix_format /= PM_index)then
        write(use_unit,*) 'Error: output_eigenvalue_eigenvector_p1 supports only packed matrix format index'
        return
     end if


     !find the node where the eigenvector is stored
     allocate(myid_tmp(n_k_points))
     
     do i_k_point = 1, n_k_points
        if (use_scalapack) then
        
          do i=0, n_scalapack_tasks-1
		  if (i_k_point == i*n_k_points/n_scalapack_tasks + 1) then
		     myid_tmp(i_k_point) = i
                 exit
		  end if
           end do
       
        else
        
           myid_tmp(i_k_point) = MOD(i_k_point, n_tasks)

       end if
    end do


    i_k = 0
    do i_k_point = 1, n_k_points
       
        		
       if (myid .eq. myid_tmp(i_k_point)) then
          i_k = i_k + 1
           
          !file name with the index of k_point appended at the end
          if (i_k_point < 10) then
             write(c_num,'("00",I1)') i_k_point
          else if (i_k_point < 100) then
             write(c_num,'("0",I2)') i_k_point
          else if (i_k_point < 1000) then
             write(c_num,'(I3)') i_k_point
          else
             write(use_unit,*)"Error: too many files to output!"
             stop
          end if
          write(filename_new,'(A)') trim(adjustl(filename))//c_num

          !each line (except line 1 and 2) contains: 1 real number and n_basis complex numbers
          IRECL = (2 + n_basis*4)*LEN
          fileID = 50 + i_k_point
          open (UNIT=fileID, FILE=filename_new, STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=IRECL)
           
          !output the two heading lines
          iline = 1
          write(fileID, REC=iline) n_basis, n_states, n_spin, n_k_points, n_atoms
          iline = iline + 1
          write(fileID, REC=iline) basis_atom(1:n_basis)

          !output eigenvector(:, :, :, i_k_point) into the file
          do i_spin = 1, n_spin, 1
             do i_state = 1, n_states, 1

                iline = iline + 1
                if (real_eigenvectors) then
                   write(fileID, REC=iline) (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential)*hartree, &
                                     dcmplx(KS_eigenvector(1:n_basis, i_state, i_spin, i_k))
                else
                   write(fileID, REC=iline) (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential)*hartree, &
                                     KS_eigenvector_complex(1:n_basis, i_state, i_spin, i_k)
                end if
                                                                                                           
             end do
          end do

          close(fileID)

        end if
     end do


  end if

  deallocate(myid_tmp)
   
end subroutine output_eigenvalue_eigenvector_p1





subroutine output_overlap_matrix_p1 ( overlap_matrix, filename )

!  PURPOSE
!  The subroutine outputs overlap matrix for periodic systems.
!  See below for the detailed information about the output file.
!******************************************************************************************************

!The content of the output unformatted file:
!line 1:  n_basis, n_spin, n_k_points
!line 2 to (n_basis * n_spin * n_k_points + 1):  
!do i_k_point = 1, n_k_points
!   do i_spin = 1, n_spin, 1
!      do i_basis = 1, n_basis, 1
!         write: overlap_matrix(1:n_basis, i_basis, i_spin, i_k_point)
!      end do
!   end do
!end do
!
!*** Attention****
!Only the upper part of the stored overlap_matrix contains useful data
!The elements in the lower part:
!overlap_matrix(i_basis_1, i_basis_2, i_spin, i_k_point) = 0 when i_basis_1 > i_basis_2
!
!Use the Hermitian property of the overlap matrix to restore the correct values for them:
!overlap_matrix(i_basis_1, i_basis_2, i_spin, i_k_point) = dconjg(overlap_matrix(i_basis_1, i_basis_2, i_spin, i_k_point))

!******************************************************************************************************
!
!  USES

  use dimensions
  use constants
  use mpi_tasks
  use localorb_io
  use basis
  use geometry
  use species_data
  use runtime_choices
  use pbc_lists
  use synchronize_mpi
  use density_matrix_evaluation
  use scalapack_wrapper
  implicit none

!  ARGUMENTS

  real*8, dimension( n_hamiltonian_matrix_size )                :: overlap_matrix
  character*40 :: filename


!  INPUTS
!   o overlap_matrix -- overlap matrix
!   o filename -- file where results are printed
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
!  SOURCE
!

  ! local variables
  complex*16, dimension(:), allocatable :: tmp_ovlp_mat
  character*100 :: info_str

  ! counters
  integer :: i_spin, i_state
  integer :: i_basis_1, i_basis_2
  integer :: i_index
  integer :: i_k, i_k_point, i_size, i_cell
  integer :: IRECL, iline
  integer :: LEN, itmp

  INQUIRE (IOLENGTH=LEN) itmp



  ! begin work

  write(info_str,'(2X,A)') &
       ' '
  call localorb_info(info_str)
  write(info_str,'(2X,A,A)') &
       'Output overlap matrix into unformatted file: ', filename
  call localorb_info(info_str)


  if(n_periodic == 0 .and. packed_matrix_format==PM_none)then
     write(use_unit,*) 'Error: output_overlap_matrix_p1 only supports periodic systems.'
     return

  else

     if(packed_matrix_format /= PM_index)then
        write(use_unit,*) 'Error: output_overlap_matrix_p1 supports only packed matrix format index'
        return
     end if


     if (myid.eq.0) then

        !*** All the elements of overlap matrix are stored on the node 0.
        !each line (except line 1) contains: n_basis complex numbers
        IRECL = n_basis*4*LEN

        open (UNIT=50, FILE=filename, STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=IRECL)
        iline = 1
        write(50, REC=iline) n_basis, n_spin, n_k_points
        
        allocate(tmp_ovlp_mat(n_basis))

        do i_k_point = 1, n_k_points
           do i_spin = 1, n_spin, 1
              do i_basis_2 = 1, n_basis
                 tmp_ovlp_mat = (0.0d0, 0.0d0)
                 
                 !use tmp_ovlp_mat to store overlap_matrix(:, i_basis_2, i_spin, i_k_point)
                 !but not that overlap_matrix(i_basis_2+1:n_basis, i_basis_2, i_spin, i_k_point) are all forced into zero

                 do i_cell = 1,n_cells_in_hamiltonian-1
                    if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                       i_index = index_hamiltonian(1,i_cell, i_basis_2)-1
                       do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                          i_index = i_index + 1
                          i_basis_1 =  column_index_hamiltonian(i_index)

                          tmp_ovlp_mat(i_basis_1) = tmp_ovlp_mat(i_basis_1) + &
                          k_phase(i_cell,i_k_point) * overlap_matrix(i_index)

                       end do
                    end if
                 end do

                 !output overlap_matrix(1:n_basis, i_basis, i_spin, i_k_point)
                 iline = iline + 1
                 write(50, REC=iline) tmp_ovlp_mat(1:n_basis)

              end do
           end do
        end do
        
        close(50)
        deallocate(tmp_ovlp_mat)
     
     end if
  end if

end subroutine output_overlap_matrix_p1
