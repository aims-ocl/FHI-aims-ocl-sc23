!****h* FHI-aims/hdf5_output
!  NAME
!    hdf5_output
!  SYNOPSIS


module hdf5_output
!  PURPOSE
!    Collection of HDF5 output routines used during band structure calculations.
!  AUTHOR
!    Lydia Nemec, FHI-aims team, Chair for Theoretical Chemistry, Technical University Munich
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  NOTE
!    Don't let the name of this module fool you; other HDF5 output routines are
!    scattered throughout the code.
!  USES
   implicit none
!  SOURCE

! By default all members below are private
   private

! public members
   public :: output_overlap_scalapack
   public :: output_complex_overlap_matrix
   public :: output_complex_eigenvector
   public :: output_complex_EV_scalapack
   public :: output_complex_KS_OCC_scalapack
   public :: output_complex_hamiltonian_matrix

   contains

!****s* FHI-aims/output_complex_overlap_matrix
!  NAME
!   output_complex_overlapmatrix
!  SYNOPSIS
   subroutine output_complex_overlap_matrix &
      ( overlap_matrix_w_complex, i_band, i_k_point, &
        i_x, i_y, i_z )

!  PURPOSE
!  For a given k-point writes the complex overlap matrix
!  to a given output file. Intention is that this subroutine be called only from MPI task
!  zero, the usual output thread.
!
!  This subroutine works with the LAPACK band structure calculation. For every band segment
!  and every k-point in the band structure calculation the overlap matrix is written in
!  either a text file or alternatively a HDF5.
!  There are to ways to activate the overlap matrix output
!  1.) specify in the control.in file:
!      output overlap_matrix
!  2.) if a keyword leading with
!      unfold
!      is given in the control.in file
!  In the case of a LAPACK band structure calculation the default output is in a text file.
!  However, by specifying in the control.in file:
!  unfold format hdf5 (for hdf5 output)
!  or
!  unfold format ASCII (for text output)
!
!  In the cse of a ScaLAPACK band structure calculation only HDF5  output is supported.
!
!  USES
   use dimensions , only: n_basis, n_k_points
   use HDF5
   use hdf5_tools, only: out_ovlp_scalapack
   use generate_aims_uuid, only: write_aims_uuid
   use localorb_io, only: localorb_info
   use mpi_tasks, only: myid, n_tasks, check_allocation
   use runtime_choices, only : flag_hdf5_unfold, flag_ASCII_unfold
   implicit none

!  ARGUMENTS
   complex*16, dimension(n_basis*(n_basis+1)/2), intent(in) :: overlap_matrix_w_complex
   integer, intent(in) :: i_band, i_k_point
   real*8, intent(in)  :: i_x,i_y,i_z

!  INPUTS
!   o i_band             -- number of present requested output band
!   o i_k_point          -- number of k-point in present band
!   o i_x, i_y, i_z      -- real(!), relative coordinates of present k-point in units
!                           of the reciprocal lattice vector of the present structure
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE

!  variables needed for HDF5 setup
   integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
   integer(HID_T)  :: plist_id  ! Property list identifier for HDF5 usage
   integer(HID_T)  :: dset_id   ! Dataset identifier
   integer(HID_T)  :: dspace_id ! Dataspace identifier
   integer(HID_T)  :: attr_id   ! Attribute identifier
   integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
   integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier

   integer(HSIZE_T), dimension(4) :: dim_ovl ! Dimensions of the overlap matrix
   integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension

   integer         :: drank = 4 ! Dataset rank
   integer         :: arank = 1 ! Attribure rank
   integer(SIZE_T) :: attrlen   ! Length of the attribute string
   character*50    :: file_name ! File name
   character*150, dimension(3) ::  attr_data  ! Attribute data
   character*10    :: num_char
   integer         :: error     ! Error flag

!  local variables    
   integer,dimension(:,:),allocatable :: index_b
   complex*16  :: overlap_matrix_entry
   character*150 :: info_str
!  counters
   integer :: i_row, i_col
   integer :: i_entry, i_index

!-----------------------------------------------------------------------!
!  begin work
!-----------------------------------------------------------------------!
   if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then

!-----------------------------------------------------------------------!
! Creating the mapping array, only needed for ASCII output
!-----------------------------------------------------------------------!
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

!-----------------------------------------------------------------------!
! Create generic filename for overlap matrix output. Tedious but general.
!-----------------------------------------------------------------------!
      write(file_name, '(A)') 'KS_overlap_matrix.'

      ! removal of trailing characters.
      write(file_name, '(A)') trim(file_name)

      ! encode band number
      write(num_char,'(I10)') i_band
      write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'

      ! encode k-point number
      write(num_char,'(I10)') i_k_point
      if (flag_hdf5_unfold) then
         write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'
      else 
         write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.out'
      end if

!-----------------------------------------------------------------------!
! Start with the writing process
!-----------------------------------------------------------------------!
      if (flag_hdf5_unfold) then
         ! Define the dimensions of the overlap matrix for HDF5 output
         ! dim_ovl has rank 4
         dim_ovl(1:4) = (/ n_basis, n_basis, 1, 2 /)
         ! Length of the attribute string
         attrlen = 150

!-----------------------------------------------------------------------!
! Opening H5 data file
!-----------------------------------------------------------------------!
! Initialize HDF5 FORTRAN interface
        call h5open_f(error)
! create the HDF5 file
        call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
! Create the dataspace.
        call h5screate_simple_f(drank, dim_ovl, dspace_id, error)
! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'cplx_ovlp_mtrx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        
!-----------------------------------------------------------------------!
! Create scalar data space for the attribute.
!-----------------------------------------------------------------------!
        call h5screate_simple_f(arank, adims, aspace_id, error)
! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        call h5tset_size_f(atype_id, attrlen, error)
! Create dataset attribute.
        call write_aims_uuid(info_str)
        attr_data(1) = info_str

        write(info_str,'(A,I8)') &
         '# Complex overlap matrix for band number ', i_band
        attr_data(2) = info_str

        write(info_str,'(A,I6,A,3(F11.8,1X))') &
         '# k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', i_x, i_y, i_z
        attr_data(3) = info_str
      
        call h5acreate_f(dset_id, "COM attribute", atype_id, aspace_id, attr_id, error)
! Write the attribute data.  
        call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
! Close the attribute.
        call h5aclose_f(attr_id, error)
!-----------------------------------------------------------------------!

! End access to the dataset and release resources used by it.
        call h5dclose_f(dset_id, error)

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
       
            call out_ovlp_scalapack( &
                 overlap_matrix_entry, &
                 file_id, &
                 plist_id, &
                 i_col, i_row, &
                 1, &
                 'cplx_ovlp_mtrx')

          end do ! end loop i_row
        end do ! end loop i_col
!-----------------------------------------------------------------------!
! Closing H5 data file
!-----------------------------------------------------------------------!
! Terminate access to the data space.
        call h5sclose_f(dspace_id, error)
! Terminate access to the file.
        call h5fclose_f(file_id, error)
! Close FORTRAN interface.!
        call h5close_f(error)

      else
!-----------------------------------------------------------------------!
! File header and initial per-state information first
!-----------------------------------------------------------------------!
        open (15, FILE=file_name)

        write(15,'(A,I8,A,I8,A,3(F12.8,1X))') &
         '# Complex overlap matrix for band number ', i_band, &
         ', k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', &
         i_x, i_y, i_z
        call write_aims_uuid(info_str)
        write(15,'(A)') info_str

        write(15,'(A)') '#'

!-----------------------------------------------------------------------!
!    find the right entries for every line in the n_basis x n_basis overlap matrix
!-----------------------------------------------------------------------!

        do i_row = 1, n_basis
          do i_col = 1, n_basis
            overlap_matrix_entry = (0.0,0.0)
            i_entry = 0
            if (i_row < i_col) then
              i_entry = index_b(i_row, i_col)
              overlap_matrix_entry = overlap_matrix_w_complex(i_entry)
              write(15,'(3X,F13.8,1X,F13.8)',ADVANCE='NO') &
                        REAL( overlap_matrix_entry ), &
                        AIMAG( overlap_matrix_entry )

            else
              i_entry = index_b(i_col, i_row)
              overlap_matrix_entry = overlap_matrix_w_complex(i_entry)
              write(15,'(3X,F13.8,1X,F13.8)',ADVANCE='NO') &
                         REAL( overlap_matrix_entry ), &
                         AIMAG(CONJG( overlap_matrix_entry ))
            end if
          end do ! end loop i_col
          write(15,'()') ! end line
        end do ! end loop i_row
        close(15)
      end if ! flag HDF output format

    write (info_str,'(2X,A)') &
        "| Writing complex overlap matrix in file:"
    call localorb_info(info_str)
    write (info_str,'(2X,A,A)') &
        "| ", file_name
    call localorb_info(info_str)
    
    deallocate(index_b)
    end if ! my_id

  end subroutine output_complex_overlap_matrix
!----------------------------------------------------------------------
!******

!-----------------------------------------------------------------------------------
!****s* output_overlap_scalapack
!  NAME
!    output_overlap_scalapack
!  SYNOPSIS
  subroutine output_overlap_scalapack( i_band, i_k_point )
!  PURPOSE
!  For a given k-point writes the complex overlap matrix into a HDF5 file.
!  This routine is based on parallel writting of all nodes.
!
!  This subroutine works with the SCALAPACK band structure calculation. For every 
!  band segment and every k-point in the band structure calculation the overlap
!  matrix is written in a HDF5.
!  There are to ways to activate the overlap matrix output
!  1.) specify in the control.in file:
!      output overlap_matrix
!  2.) if a keyword leading with
!      unfold
!      is given in the control.in file
!  In the case of a SCALAPACK band structure calculation only HDF5 output
!  is implemented.
!
!  USES
    use dimensions , only: n_basis, n_k_points, n_hamiltonian_matrix_size
    use HDF5
    use hdf5_tools, only: open_hdf5, close_hdf5, out_ovlp_scalapack
    use generate_aims_uuid, only: write_aims_uuid
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: aims_stop, mpi_comm_global
    use pbc_lists
    use runtime_choices, only: flag_hdf5_unfold, flag_ASCII_unfold, &
        use_local_index, packed_matrix_format, n_points_in_band, band_begin, &
        band_end, PM_index
    use scalapack_wrapper, only: my_k_point, l_col, l_row, ovlp_complex
   implicit none

!  ARGUMENTS
    integer, intent(in) :: i_band, i_k_point

!  variables needed for HDF5 setup
    integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
    integer(HID_T)  :: plist_id  ! Property list identifier for HDF5 usage
    integer(HID_T)  :: dset_id   ! Dataset identifier
    integer(HID_T)  :: dspace_id ! Dataspace identifier
    integer(HID_T)  :: attr_id   ! Attribute identifier
    integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
    integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier
    
    integer(HSIZE_T), dimension(4) :: dim_ovl ! Dimensions of the overlap matrix
    integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension

    
    integer         :: drank = 4 ! Dataset rank
    integer         :: arank = 1 ! Attribure rank   
    integer(SIZE_T) :: attrlen   ! Length of the attribute string
    character*50    :: file_name ! File name
    character*150, dimension(3) ::  attr_data  ! Attribute data
    character*10    :: num_char
    integer         :: error     ! Error flag
    
! MPI Flag
    integer :: info

    character*150 :: info_str

    
    integer:: i_k, i_cell, i_col, i_row, lr, lc, idx
    integer:: i_k_point_band
    real*8:: i_x, i_y, i_z


    
! Define the dimensions of the overlap matrix for HDF5 output
! dim_ovl has rank 4 
    dim_ovl(1:4) = (/ n_basis, n_basis, 1, 2 /)
! Length of the attribute string
    attrlen = 150
    
    
    if(use_local_index) call aims_stop("INTERNAL ERROR: output_overlap_scalapack + use_local_index")
    
    do i_k = 1, n_k_points
      i_k_point_band = i_k_point+i_k-1
      
      if(i_k_point_band > n_points_in_band(i_band)) cycle

! Calculate reciprocal coordinates of current k-point
      i_x =  band_begin(i_band,1) +  real(i_k-1)/real(n_points_in_band(i_band)-1) &
             *( band_end(i_band,1) -  band_begin(i_band,1))
      i_y =  band_begin(i_band,2) +  real(i_k-1)/real(n_points_in_band(i_band)-1) &
             *( band_end(i_band,2) -  band_begin(i_band,2))
      i_z =  band_begin(i_band,3) +  real(i_k-1)/real(n_points_in_band(i_band)-1) &
             *( band_end(i_band,3) -  band_begin(i_band,3))
             
!-----------------------------------------------------------------------!
! Create generic filename for overlap matrix output. Tedious but general.
!-----------------------------------------------------------------------!
      write(file_name, '(A)') 'KS_overlap_matrix.'

! removal of trailing characters.
      write(file_name, '(A)') trim(file_name)

! encode band number
      write(num_char,'(I10)') i_band
      write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'

! encode k-point number
      write(num_char,'(I10)') i_k_point_band
      write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'

      write (info_str,'(2X,2A)') '| Writing complex overlap matrix: File ', file_name
      call localorb_info ( info_str ) 
!-----------------------------------------------------------------------!

! Create HDF5 file.
      call open_hdf5(file_name, file_id, plist_id)
      call mpi_barrier(mpi_comm_global,info)
       
! Create the dataspace.
      call h5screate_simple_f(drank, dim_ovl, dspace_id, error)
      call h5dcreate_f(file_id, 'cplx_ovlp_mtrx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
              
!-----------------------------------------------------------------------!
! Create scalar data space for the attribute.
!-----------------------------------------------------------------------!
      call h5screate_simple_f(arank, adims, aspace_id, error)
! Create datatype for the attribute.
      call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
! Create dataset attribute.
      call write_aims_uuid(info_str)
      attr_data(1) = info_str

      write(info_str,'(A,I8)') &
         '# Complex overlap matrix for band number ', i_band
      attr_data(2) = info_str

      write(info_str,'(A,I6,A,3(F11.8,1X))') &
         '# k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', i_x, i_y, i_z
      attr_data(3) = info_str
      
      call h5acreate_f(dset_id, "COM attribute", atype_id, aspace_id, attr_id, error)
! Write the attribute data.  
      call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
! Close the attribute.
      call h5aclose_f(attr_id, error)
!-----------------------------------------------------------------------!
      call h5dclose_f(dset_id, error)

      select case(packed_matrix_format)

        case(PM_index)
          if(i_k==my_k_point)then
      
            do i_cell = 1, n_cells_in_hamiltonian-1
              do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                  do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                    i_row = column_index_hamiltonian(idx)
                    lr = l_row(i_row) ! local row number
                    if(lr==0) cycle   ! skip if not local
                     
                    call out_ovlp_scalapack( &
                       ovlp_complex(lr,lc), &
                       file_id, &
                       plist_id, &
                       i_col, i_row, &
                       1, &
                       'cplx_ovlp_mtrx')

                    call out_ovlp_scalapack( &
                       CONJG(ovlp_complex(lr,lc)), &
                       file_id, &
                       plist_id, &
                       i_row, i_col, &
                       1, &
                       'cplx_ovlp_mtrx')
                  end do  ! idx
                end if  ! index_hamiltonian(1,i_cell,i_col) > 0
              end do  ! i_coll
            end do  ! i_cell
          end if  ! i_k==my_k_point

      case default

        write(use_unit,*) 'Error: construct_overlap_scalapack does not support non-packed matrices.'
        call aims_stop

      end select ! packed_matrix_format

      call mpi_barrier(mpi_comm_global,info)

!-----------------------------------------------------------------------!
! Closing H5 data file
!-----------------------------------------------------------------------!
! Terminate access to the data space.
      call h5sclose_f(dspace_id, error)

      call close_hdf5(file_id,plist_id)

      write (info_str,'(2X,2A)') '| ... finished writing complex overlap matrix: File ', file_name
      call localorb_info ( info_str ) 
    
    end do  ! i_k

  end subroutine output_overlap_scalapack
  !******

!****s* FHI-aims/output_complex_eigenvector
!  NAME
!   output_complex_eigenvector
!  SYNOPSIS

  subroutine output_complex_eigenvector &
      ( KS_eigenvector_tmp, KS_eigenvalue, chemical_potential, occ_numbers, & 
        i_band, i_k_point, i_spin, i_x, i_y, i_z &
      )

!  PURPOSE
!
!  This subroutine works with the LAPACK band structure calculation.
!  For every band segment, spin channel, and k-point in the
!  band structure calculation the complex Kohn-Sham eigenvector
!  c_il(k) (i=basis, l=state) are written in either a text file
!  or alternatively a HDF5. Intention is that this subroutine be
!  called only from MPI task zero, the usual output thread.
!  There are to ways to activate the eigenvector output
!  1.) specify in the control.in file:
!      output eigenvectors
!  2.) if a keyword leading with
!      unfold
!      is given in the control.in file
!  In the case of a LAPACK band structure calculation the default output is in a text file.
!  However, by specifying in the control.in file:
!  unfold format hdf5 (for hdf5 output)
!  or
!  unfold format ASCII (for text output)
!  the output format can be specified.
!
!  USES
    use basis, only : basis_fn, basis_atom, basisfn_type, basisfn_n, basis_m, &
        basis_l
    use constants, only: hartree
    use dimensions, only : n_basis, n_states, n_k_points, n_spin
    use HDF5
    use generate_aims_uuid, only: write_aims_uuid
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: myid, check_allocation
    use runtime_choices, only : flag_hdf5_unfold, flag_ASCII_unfold
    implicit none

!  ARGUMENTS

    complex*16, dimension(n_basis, n_states), intent(in) :: KS_eigenvector_tmp
    real*8, dimension(n_states), intent(in)              :: KS_eigenvalue
    real*8, intent(in)                                   :: chemical_potential
    real*8, dimension(n_states), intent(in)              :: occ_numbers
    integer, intent(in) :: i_band, i_k_point, i_spin
    real*8, intent(in)  :: i_x,i_y,i_z

!  INPUTS
!   o KS_eigenvector_tmp -- Kohn-Sham eigenvectors, only for one k-point
!   o KS_eigenvalue      -- Kohn-Sham eigenvalues for present k-point
!   o chemical_potential -- Fermi level
!   o occ_numbers        -- occupation weights of eigenstates for present k-point
!   o i_band             -- number of present requested output band
!   o i_k_point          -- number of k-point in present band
!   o i_spin             -- number of present spin
!   o i_x, i_y, i_z      -- real(!), relative coordinates of present k-point in units
!                           of the reciprocal lattice vector of the present structure
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE


!  local variables
    character l_char
    character*150 :: info_str
    real*8, dimension(:),allocatable :: KS_state_eV
    real*8, dimension(:,:,:),allocatable :: EV_entry
    integer, dimension(:,:),allocatable :: basis_list
!  counters
    integer i_basis, i_state, i_fn, info

!  functions
    character l_to_str

!  variables needed for HDF5 setup
    integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
    integer(HID_T)  :: plist_id  ! Property list identifier for HDF5 usage
    integer(HID_T)  :: dset_id   ! Dataset identifier
    integer(HID_T)  :: dspace_id ! Dataspace identifier
    integer(HID_T)  :: attr_id   ! Attribute identifier
    integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
    integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier
    integer(HID_T)  :: memspace  ! Memory space identifier

    integer(HSIZE_T), dimension(1) :: dim_states ! Dimensions of the eigenvalue array
    integer(HSIZE_T), dimension(1) :: dim_basis  ! Dimensions of the basis function ID array
    integer(HSIZE_T), dimension(3) :: dim_Evec   ! Dimensions of the eigenvector array    
    integer(HSIZE_T), dimension(3) :: dim_entry  ! Dimensions of the eigenvector entry  
    integer(HSIZE_T), dimension(2) :: dim_B_list ! Dimensions of the entry containing the basis function angular momentum
    integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension
    integer(HSIZE_T), dimension(3) :: offset
    
    integer         :: drank = 1 ! Dataset rank
    integer         :: dEVrank = 3 ! Dataset rank
    integer         :: arank = 1 ! Attribure rank
    integer(SIZE_T) :: attrlen   ! Length of the attribute string
    character*150, dimension(3) ::  attr_data  ! Attribute data
    character*10    :: num_char
    character*50    :: file_name
    integer         :: error     ! Error flag

!-----------------------------------------------------------------------!
!  begin work
!-----------------------------------------------------------------------!


      !-----------------------------------------------------------------!
      ! Generate the output file name
      !-----------------------------------------------------------------!
      ! spin naming first
      if (n_spin.eq.1) then
        write(file_name, '(A)') 'KS_eigenvectors.'
      else if (i_spin.eq.1) then
        write(file_name, '(A)') 'KS_eigenvectors_up.'
      else
        write(file_name, '(A)') 'KS_eigenvectors_dn.'
      end if
      !
      ! encode band number
      write(num_char,'(I10)') i_band
      write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'
      ! encode k-point number
      write(num_char,'(I10)') i_k_point
      if (flag_hdf5_unfold) then
        write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'
      else
         write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.out'
      end if
      !
      ! removal of trailing characters.
      write(file_name, '(A)') trim(file_name)
      !-------------------------------------------------------------------!

    if (myid.eq.0) then

      write (info_str,'(2X,A)') &
        "| Writing complex Kohn-Sham eigenvectors in file:"
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(2X,A,A)') &
        "| ", file_name
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'()') ! end line
      call localorb_info(info_str,use_unit,'(A)')
        
      if (flag_hdf5_unfold) then
          
          ! Define the dimensions of the basis function ID for HDF5 output
          ! dim_basis has rank 1
          dim_basis(1) = n_basis
          
          ! Define the dimensions of the eigenstates and occ_numbers for HDF5 output
          ! dim_states has rank 1
          dim_states(1) = n_states
          
        !-----------------------------------------------------------------!
        ! Opening H5 data file
        !-----------------------------------------------------------------!
        ! Initialize HDF5 FORTRAN interface
        call h5open_f(error)
        ! create the HDF5 file
        call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error)

        ! Generate array with all eigenvalues in eV for HDF5 writing
        allocate(KS_state_eV(n_states), stat=info)
        call check_allocation(info, 'KS_state_eV')

        do i_state= 1,n_states,1
          KS_state_eV(i_state) = (KS_eigenvalue(i_state)-chemical_potential)*hartree
        end do
        
        ! Write eigenstates to dataset
        call h5screate_simple_f(drank, dim_states, dspace_id, error)
        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'KS_EigenStates', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
              
        ! Write KS_state_eV to HDF5 dataset and allocate the array as it is no longer needed
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, KS_state_eV, dim_states, error)
        deallocate(KS_state_eV)
      
        ! End access to the dataset and release resources used by it.
        call h5dclose_f(dset_id, error)
        ! Terminate access to the data space.
        call h5sclose_f(dspace_id, error)
        !-----------------------------------------------------------------!

        !-----------------------------------------------------------------!
        ! Write occupation numbers to dataset
        ! create the dataspace.
        call h5screate_simple_f(drank, dim_states, dspace_id, error)

        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'KS_occ_numbers', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

        ! Write occupation numbers to HDF5 dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, occ_numbers, dim_states, error)                      

        ! End access to the dataset and release resources used by it.
        call h5dclose_f(dset_id, error)
        ! Terminate access to the data space.
        call h5sclose_f(dspace_id, error)
        !-----------------------------------------------------------------!

        !-----------------------------------------------------------------!
        ! Write KS eigenvector to dataset
        dim_Evec(1:3) = (/ n_states, n_basis, 2 /)
        
        ! create the dataspace.
        call h5screate_simple_f(dEVrank, dim_Evec, dspace_id, error)
        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'cplx_eigenvect', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

        !-----------------------------------------------------------------!
        ! Create scalar data space for the attribute.
        !-----------------------------------------------------------------!
        ! Length of the attribute string
        attrlen = 150
        call h5screate_simple_f(arank, adims, aspace_id, error)
        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        call h5tset_size_f(atype_id, attrlen, error)
        ! Create dataset attribute.
        call write_aims_uuid(info_str)
        attr_data(1) = info_str
        write(info_str,'(A,I8)') &
         '# Complex Kohn-Sham eigenvectors for band number ', i_band
        attr_data(2) = info_str

        write(info_str,'(A,I6,A,3(F11.8,1X))') &
         '# k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', i_x, i_y, i_z
        attr_data(3) = info_str
      
        call h5acreate_f(dset_id, "COM attribute", atype_id, aspace_id, attr_id, error)
        ! Write the attribute data.
        call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
        !-----------------------------------------------------------------!
        dim_entry = (/1,1,2/)
        allocate(EV_entry(1,1,2))

        do i_basis = 1, n_basis, 1
          do i_state = 1, n_states, 1
            
            !!dim_Evec(1:3) = (/ n_states, n_basis, 2 /)
            EV_entry(1,1,1)=REAL(KS_eigenvector_tmp(i_basis,i_state))
            EV_entry(1,1,2)=AIMAG(KS_eigenvector_tmp(i_basis,i_state))
            
            offset(1:3) = (/ i_state-1, i_basis-1, 0 /)
            call h5screate_simple_f(dEVrank, dim_entry, memspace, error)
            call h5dget_space_f(dset_id, dspace_id, error)
            call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, &
                                        dim_entry, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, EV_entry, dim_entry, &
                            error, file_space_id = dspace_id, &
                            mem_space_id = memspace)
            call h5sclose_f(memspace, error)
          end do
        end do
        deallocate(EV_entry)
        ! End access to the eigenvector dataset and release resources used by it.
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        !-----------------------------------------------------------------!
        ! Write Basis function informations to dataset
        !-----------------------------------------------------------------!
        ! create the dataspace.
        allocate(basis_list(n_basis,4), stat=info)
        call check_allocation(info, 'basis_list')        
        
        do i_basis = 1, n_basis, 1
          i_fn = basis_fn(i_basis)
          
          basis_list(i_basis,1:4) = (/ basis_atom(i_basis), &
                                       basisfn_n(i_fn),     &
                                       basis_l(i_basis),    &
                                       basis_m(i_basis)     /)
        end do
        
        dim_B_list(1:2) = (/ n_basis, 4 /)
        
        call h5screate_simple_f(2, dim_B_list, dspace_id, error)
        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'Basis_funct_ID', H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        ! Write occupation numbers to HDF5 dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, basis_list, dim_B_list, error)
        ! End access to the dataset and release resources used by it.
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        
        deallocate(basis_list)

        !-----------------------------------------------------------------!
        ! Closing H5 data file
        !-----------------------------------------------------------------!

        ! Terminate access to the file.
        call h5fclose_f(file_id, error)
        ! Close FORTRAN interface!
        call h5close_f(error)
      else
        !-----------------------------------------------------------------!
        ! Start writing data to ASCII file
        !-----------------------------------------------------------------!
        open (50, FILE=file_name)

        ! File header and initial per-state information first
        call write_aims_uuid(info_str)
        write(50,'(A)') info_str

        write(50,'(A,I8,A,I8,A,3(F12.8,1X))') &
            '# Complex Kohn-Sham eigenvectors for band number ', i_band, &
            ', k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', &
            i_x, i_y, i_z
        write(50,'(A)') '#'
        write(50,'(A,6X,A,23X,A)') '#','Basis function ID','State number'
        write(50,'(A,4X,A,2X)',ADVANCE='NO') '#','/                   \ '

        do i_state = 1, n_states, 1
          write(50,'(I30)',ADVANCE='NO') i_state
        enddo
        write(50,'()') ! end line

        write(50,'(A,3X,A,1X)',ADVANCE='NO') '#','/                     \ '

        do i_state=1,n_states,1
          write(50,'(10X,F20.8)',ADVANCE='NO') (KS_eigenvalue(i_state)-chemical_potential)*hartree
        end do
        write(50,'()') ! end line

        write(50,'(A,2X,A)',ADVANCE='NO')    '#','/                       \ '

        do i_state=1,n_states,1
        write(50,'(10X,F20.8)',ADVANCE='NO') occ_numbers(i_state)
        end do
        write(50,'()') ! end line

        write(50,'(A)',ADVANCE='NO') '# No.  atom   type    n l   m'

        do i_state=1,n_states,1
          write(50,'(5X,A11,3X,A11)',ADVANCE='NO') 'Re[c_il(k)]', 'Im[c_il(k)]'
        end do
        write(50,'()') ! end line
        !-----------------------------------------------------------------!
        ! then follow the actual eigenvectors (real and imaginary parts)
        ! for each basis function
        !-----------------------------------------------------------------!
        do i_basis = 1, n_basis, 1

          i_fn = basis_fn(i_basis)
          l_char = l_to_str(basis_l(i_basis))

          write(50, &
              '(I5,1X,I5,1X,A8,I3,1X,A1,1X,I3)',ADVANCE='NO') &
               i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
               basisfn_n(i_fn), l_char, basis_m(i_basis)

          do i_state = 1, n_states, 1
            write(50,'(3X,F13.8,1X,F13.8)',ADVANCE='NO') &
                REAL(KS_eigenvector_tmp(i_basis,i_state)), &
                AIMAG(KS_eigenvector_tmp(i_basis,i_state))
          end do
          write(50,'()') ! end line

        end do ! i_basis

        close(50)

      end if  ! flag HDF output format

    end if ! my_id
    
    end subroutine output_complex_eigenvector
!----------------------------------------------------------------------
!******

!****s* FHI-aims/output_complex_EV_scalapack
!  NAME
!   output_complex_EV_scalapack
!  SYNOPSIS

  subroutine output_complex_EV_scalapack ( i_band, i_k_point, i_spin)

!  PURPOSE
!
!  This subroutine works with the SCALAPACK band structure calculation.
!  For every band segment, spin channel, and k-point in the
!  band structure calculation the complex Kohn-Sham eigenvector
!  c_il(k) (i=basis, l=state) are written in a HDF5. The Eigenvalue, occupation numbers
!  and Basis set information will be written from an extra subroutine.
!
!  Intention is that this subroutine be called only from MPI task zero,
!  the usual output thread.
!
!  There are to ways to activate the eigenvector output
!  1.) specify in the control.in file:
!      output eigenvectors
!  2.) if a keyword leading with
!      unfold
!      is given in the control.in file
!  In the case of a SCALAPACK band structure calculation the default output is in a HDF5 file.
!  Output in a textfile is not supported for SCALAPACK EV output.
!
!  USES

    use basis, only : basis_fn, basis_atom, basisfn_type, basisfn_n, basis_m, &
        basis_l
    use constants, only: hartree
    use dimensions, only : n_basis, n_states, n_k_points, n_spin
    use generate_aims_uuid, only: write_aims_uuid
    use HDF5
    use hdf5_tools, only: open_hdf5, close_hdf5
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: mpi_comm_global, aims_stop
    use runtime_choices, only: flag_hdf5_unfold, flag_ASCII_unfold, &
        use_local_index, packed_matrix_format, n_points_in_band, band_begin, &
        band_end, PM_index
    use scalapack_wrapper, only: my_k_point, eigenvec_complex, l_row, l_col
    implicit none

!  ARGUMENTS

    integer, intent(in) :: i_band, i_k_point, i_spin

!  INPUTS
!   o i_band             -- number of present requested output band
!   o i_k_point          -- number of k-point in present band
!   o i_spin             -- number of present spin

!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE


!  local variables
    character l_char
    character*150 :: info_str
    real*8, dimension(:,:,:),allocatable :: EV_entry
    real*8  i_x,i_y,i_z
!  counters
    integer i_basis, i_state, info, i_k_point_band, i_k

!  variables needed for HDF5 setup
    integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
    integer(HID_T)  :: plist_id  ! Property list identifier for HDF5 usage
    integer(HID_T)  :: dset_id   ! Dataset identifier
    integer(HID_T)  :: dspace_id ! Dataspace identifier
    integer(HID_T)  :: attr_id   ! Attribute identifier
    integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
    integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier
    integer(HID_T)  :: memspace  ! Memory space identifier

    integer(HSIZE_T), dimension(3) :: dim_Evec   ! Dimensions of the eigenvector array
    integer(HSIZE_T), dimension(3) :: dim_entry  ! Dimensions of the eigenvector entry
    integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension
    integer(HSIZE_T), dimension(3) :: offset

    integer         :: dEVrank = 3 ! Dataset rank
    integer         :: arank = 1 ! Attribure rank
    integer(SIZE_T) :: attrlen   ! Length of the attribute string
    character*150, dimension(3) ::  attr_data  ! Attribute data
    character*10    :: num_char
    character*50    :: file_name
    integer         :: error     ! Error flag

!-----------------------------------------------------------------------!
!  begin work
!-----------------------------------------------------------------------!
    if(use_local_index) call aims_stop("INTERNAL ERROR: output_complex_EV_scalapack + use_local_index")
    
    do i_k = 1, n_k_points
      i_k_point_band = i_k_point+i_k-1
      
      if(i_k_point_band > n_points_in_band(i_band)) cycle
      
        ! Calculate reciprocal coordinates of current k-point
        i_x =  band_begin(i_band,1) +  real(i_k_point_band-1)/real(n_points_in_band(i_band)-1) &
               *( band_end(i_band,1) -  band_begin(i_band,1))
        i_y =  band_begin(i_band,2) +  real(i_k_point_band-1)/real(n_points_in_band(i_band)-1) &
               *( band_end(i_band,2) -  band_begin(i_band,2))
        i_z =  band_begin(i_band,3) +  real(i_k_point_band-1)/real(n_points_in_band(i_band)-1) &
               *( band_end(i_band,3) -  band_begin(i_band,3))
        !-----------------------------------------------------------------!
        ! Generate the output file name
        !-----------------------------------------------------------------!
        ! spin naming first
        if (n_spin.eq.1) then
          write(file_name, '(A)') 'KS_eigenvectors.'
        else if (i_spin.eq.1) then
          write(file_name, '(A)') 'KS_eigenvectors_up.'
        else
          write(file_name, '(A)') 'KS_eigenvectors_dn.'
        end if
        !
        ! encode band number
        write(num_char,'(I10)') i_band
        write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'
        ! encode k-point number
        write(num_char,'(I10)') i_k_point_band
        if (flag_hdf5_unfold) then
          write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'
        else
           write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.out'
        end if
        !
        ! removal of trailing characters.
        write(file_name, '(A)') trim(file_name)
        
        write (info_str,'(2X,A)') &
               "| Writing complex Kohn-Sham eigenvectors in file:"
        call localorb_info(info_str,use_unit,'(A)')

        write (info_str,'(2X,A,A)') &
               "| ", file_name
        call localorb_info(info_str,use_unit,'(A)')

        write(info_str,'()') ! end line
        call localorb_info(info_str,use_unit,'(A)')
        !-------------------------------------------------------------------!
        
        ! Define the dimensions of the Eigenvector for HDF5 output
        ! dim_Evec has rank 3
         dim_Evec(1:3) = (/ n_states, n_basis, 2 /)
         dim_entry = (/1,1,2/)
         
        !-----------------------------------------------------------------!
        ! Opening H5 data file
        !-----------------------------------------------------------------!
        ! Create HDF5 file.
        call open_hdf5(file_name, file_id, plist_id)
        call mpi_barrier(mpi_comm_global,info)
        !-----------------------------------------------------------------!

        !-----------------------------------------------------------------!
        ! Write KS eigenvector to dataset
        ! Create the dataspace.
        call h5screate_simple_f(dEVrank, dim_Evec, dspace_id, error)
        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'cplx_eigenvect', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        
        !-----------------------------------------------------------------!
        ! Create scalar data space for the attribute.
        !-----------------------------------------------------------------!
        ! Length of the attribute string
        attrlen = 150
        call h5screate_simple_f(arank, adims, aspace_id, error)
        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        call h5tset_size_f(atype_id, attrlen, error)
        ! Create dataset attribute.
        call write_aims_uuid(info_str)
        attr_data(1) = info_str

        write(info_str,'(A,I8)') &
         '# Complex Kohn-Sham eigenvectors for band number ', i_band
        attr_data(2) = info_str

        write(info_str,'(A,I6,A,3(F11.8,1X))') &
         '# k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', i_x, i_y, i_z
        attr_data(3) = info_str

        call h5acreate_f(dset_id, "CEV attribute", atype_id, aspace_id, attr_id, error)
        ! Write the attribute data.
        call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
        !-----------------------------------------------------------------!
        
        select case(packed_matrix_format)

          case(PM_index)
        
            if(i_k==my_k_point) then
              allocate(EV_entry(1,1,2))
              
              call h5screate_simple_f(dEVrank, dim_entry, memspace, error)
              call h5dget_space_f(dset_id, dspace_id, error)
              
              do i_state = 1, n_states
                if(l_col(i_state)==0) cycle
                  do i_basis = 1, n_basis
                    if(l_row(i_basis)>0) then

                    call h5screate_simple_f(dEVrank, dim_entry, memspace, error)
                      call h5dget_space_f(dset_id, dspace_id, error)

                      offset(1:3) = (/ i_state-1, i_basis-1, 0 /)
                      
                      EV_entry(1,1,1)= REAL(eigenvec_complex(l_row(i_basis), l_col(i_state), i_spin))
                      EV_entry(1,1,2)=AIMAG(eigenvec_complex(l_row(i_basis), l_col(i_state), i_spin))

                      call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, &
                                                  dim_entry, error)
                      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, EV_entry, dim_entry, &
                                      error, file_space_id = dspace_id, &
                                      mem_space_id = memspace)
                      
                    end if
                  end do ! i_basis
              end do ! i_state
              
              call h5sclose_f(memspace, error)
              deallocate(EV_entry)
            end if  ! i_k==my_k_point
            
          case default
          
            call aims_stop("INTERNAL Error: output_complex_EV_scalapack does not support non-packed matrices.")

        end select ! packed_matrix_format

        call mpi_barrier(mpi_comm_global,info)
        ! End access to the eigenvector dataset and release resources used by it.
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        
        !-----------------------------------------------------------------!
        ! Closing H5 data file
        !-----------------------------------------------------------------!
        call close_hdf5(file_id,plist_id)

    end do  ! i_k

    end subroutine output_complex_EV_scalapack
!----------------------------------------------------------------------
!******

!****s* FHI-aims/output_complex_KS_OCC_scalapack
!  NAME
!   output_complex_KS_OCC_scalapack
!  SYNOPSIS

  subroutine output_complex_KS_OCC_scalapack &
      ( KS_eigenvalue, chemical_potential, occ_numbers, &
        i_band, i_k_point, i_spin)

!  PURPOSE
!
!  This subroutine works with the SCALAPACK band structure calculation.
!  The Eigenvalue, occupation numbers and Basis set information will be written
!  into the Eigenvector output HDF5 file for every band segment, spin channel, and k-point in the
!  band structure calculation. Intention is that this subroutine be called only from MPI task zero,
!  the usual output thread. This subroutine is part of the Eigenvector output called by the keyword
!  output eigenvectors or if a keyword leading with unfold is given in the control.in file.
!  In the case of a SCALAPACK band structure calculation the default output is in a HDF5 file.
!  Output in a textfile is not supported for SCALAPACK EV output.
!
!  USES
    use basis, only : basis_fn, basis_atom, basisfn_type, basisfn_n, basis_m, &
        basis_l
    use constants, only: hartree
    use dimensions, only : n_basis, n_states, n_k_points, n_spin
    use generate_aims_uuid, only: write_aims_uuid
    use HDF5
    use hdf5_tools, only: open_hdf5
    use mpi_tasks, only: myid, check_allocation
    use runtime_choices, only : flag_hdf5_unfold, flag_ASCII_unfold
    implicit none

!  ARGUMENTS

    real*8, dimension(n_states), intent(in)              :: KS_eigenvalue
    real*8, intent(in)                                   :: chemical_potential
    real*8, dimension(n_states), intent(in)              :: occ_numbers
    integer, intent(in) :: i_band, i_k_point, i_spin

!  INPUTS
!   o KS_eigenvalue      -- Kohn-Sham eigenvalues for present k-point
!   o chemical_potential -- Fermi level
!   o occ_numbers        -- occupation weights of eigenstates for present k-point
!   o i_band             -- number of present requested output band
!   o i_k_point          -- number of k-point in present band
!   o i_spin             -- number of present spin
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE


!  local variables
    character*150 :: info_str
    real*8, dimension(:),allocatable :: KS_state_eV
    integer, dimension(:,:),allocatable :: basis_list
!  counters
    integer i_basis, i_state, i_fn, info

!  variables needed for HDF5 setup
    integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
    integer(HID_T)  :: dset_id   ! Dataset identifier
    integer(HID_T)  :: dspace_id ! Dataspace identifier
    integer(HID_T)  :: attr_id   ! Attribute identifier
    integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
    integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier

    integer(HSIZE_T), dimension(1) :: dim_states ! Dimensions of the eigenvalue array
    integer(HSIZE_T), dimension(1) :: dim_basis  ! Dimensions of the basis function ID array
    integer(HSIZE_T), dimension(2) :: dim_B_list ! Dimensions of the entry containing the basis function angular momentum
    integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension

    integer         :: drank = 1 ! Dataset rank
    integer         :: arank = 1 ! Attribure rank
    integer(SIZE_T) :: attrlen   ! Length of the attribute string
    character*150, dimension(3) ::  attr_data  ! Attribute data
    character*10    :: num_char
    character*50    :: file_name
    integer         :: error     ! Error flag

!-----------------------------------------------------------------------!
!  begin work
!-----------------------------------------------------------------------!
    !-----------------------------------------------------------------!
    ! Generate the output file name
    !-----------------------------------------------------------------!
    ! spin naming first
    if (n_spin.eq.1) then
      write(file_name, '(A)') 'KS_eigenvectors.'
    else if (i_spin.eq.1) then
      write(file_name, '(A)') 'KS_eigenvectors_up.'
    else
      write(file_name, '(A)') 'KS_eigenvectors_dn.'
    end if
    !
    ! encode band number
    write(num_char,'(I10)') i_band
    write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'
    ! encode k-point number
    write(num_char,'(I10)') i_k_point
    if (flag_hdf5_unfold) then
      write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'
    else
       write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.out'
    end if
    !
    ! removal of trailing characters.
    write(file_name, '(A)') trim(file_name)
    !-------------------------------------------------------------------!

    if (myid.eq.0) then
    
      ! Define the dimensions of the basis function ID for HDF5 output
      ! dim_basis has rank 1
      dim_basis(1) = n_basis

      ! Define the dimensions of the eigenstates and occ_numbers for HDF5 output
      ! dim_states has rank 1
      dim_states(1) = n_states

      !-----------------------------------------------------------------!
      ! Opening H5 data file
      !-----------------------------------------------------------------!
      ! Initialize HDF5 FORTRAN interface
      call h5open_f(error)
      ! create the HDF5 file
      call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

      ! Generate array with all eigenvalues in eV for HDF5 writing
      allocate(KS_state_eV(n_states), stat=info)
      call check_allocation(info, 'KS_state_eV')

      do i_state= 1,n_states,1
        KS_state_eV(i_state) = (KS_eigenvalue(i_state)-chemical_potential)*hartree
      end do

      ! Write eigenstates to dataset
      call h5screate_simple_f(drank, dim_states, dspace_id, error)
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id, 'KS_EigenStates', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

      ! Write KS_state_eV to HDF5 dataset and allocate the array as it is no longer needed
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, KS_state_eV, dim_states, error)
      deallocate(KS_state_eV)

      ! End access to the dataset and release resources used by it.
      call h5dclose_f(dset_id, error)
      ! Terminate access to the data space.
      call h5sclose_f(dspace_id, error)
      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!
      ! Write occupation numbers to dataset
      ! create the dataspace.
      call h5screate_simple_f(drank, dim_states, dspace_id, error)

      ! Create the dataset with default properties.
      call h5dcreate_f(file_id, 'KS_occ_numbers', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

      ! Write occupation numbers to HDF5 dataset
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, occ_numbers, dim_states, error)

      ! End access to the dataset and release resources used by it.
      call h5dclose_f(dset_id, error)
      ! Terminate access to the data space.
      call h5sclose_f(dspace_id, error)
      !-----------------------------------------------------------------!


      !-----------------------------------------------------------------!
      ! Write Basis function informations to dataset
      !-----------------------------------------------------------------!
      ! create the dataspace.
      allocate(basis_list(n_basis,4), stat=info)
      call check_allocation(info, 'basis_list')

      do i_basis = 1, n_basis, 1
        i_fn = basis_fn(i_basis)

        basis_list(i_basis,1:4) = (/ basis_atom(i_basis), &
                                     basisfn_n(i_fn),     &
                                     basis_l(i_basis),    &
                                     basis_m(i_basis)     /)
      end do

      dim_B_list(1:2) = (/ n_basis, 4 /)

      call h5screate_simple_f(2, dim_B_list, dspace_id, error)
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id, 'Basis_funct_ID', H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
      ! Write basis set information to HDF5 dataset
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, basis_list, dim_B_list, error)
      ! End access to the dataset and release resources used by it.
      
      !-----------------------------------------------------------------!
      ! Create scalar data space for the attribute.
      !-----------------------------------------------------------------!
      ! Length of the attribute string
      attrlen = 150
      call h5screate_simple_f(arank, adims, aspace_id, error)
      ! Create datatype for the attribute.
      call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
      ! Create dataset attribute.
      call write_aims_uuid(info_str)
      attr_data(1) = info_str

      write(info_str,'(A,I8)') &
       '# Basis set information of the eigenvectors for band number', i_band
      attr_data(2) = info_str

      write(info_str,'(A)') &
       '# No.  atom, n, l, m'
      attr_data(3) = info_str
      
      call h5acreate_f(dset_id, "CEV basis-set", atype_id, aspace_id, attr_id, error)
      ! Write the attribute data.
      call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
      ! Close the attribute.
      call h5aclose_f(attr_id, error)
      !-----------------------------------------------------------------!
        
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      deallocate(basis_list)

      !-----------------------------------------------------------------!
      ! Closing H5 data file
      !-----------------------------------------------------------------!

      ! Terminate access to the file.
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface!
      call h5close_f(error)

    end if ! my_id

    end subroutine output_complex_KS_OCC_scalapack
!----------------------------------------------------------------------
!******	

!****s* FHI-aims/output_complex_hamiltonianmatrix
!  NAME
!   output_complex_hamiltonianmatrix
!  SYNOPSIS
  subroutine output_complex_hamiltonian_matrix &
      ( hamiltonian_w_complex, i_band, i_k_point, &
        i_x, i_y, i_z )

!  PURPOSE
!  For a given k-point writes the complex hamiltonian matrix
!  to a given output file. 
!  Intention is that this subroutine be called only from MPI task
!  zero, the usual output thread. This routine only works in 
!  combination with output band.

!  In the case of a LAPACK band structure calculation the default output is in a text file.
!  However, by specifying in the control.in file:
!  unfold format hdf5 (for hdf5 output)
!  or
!  unfold format ASCII (for text output)
!  the format can be changed.
!
!  USES
    use dimensions , only: n_basis, n_k_points, n_spin
    use HDF5
    use hdf5_tools, only: out_hmlt_scalapack
    use generate_aims_uuid, only: write_aims_uuid
    use localorb_io, only: localorb_info
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use runtime_choices, only : flag_hdf5_unfold, flag_ASCII_unfold
    
    implicit none

!  ARGUMENTS
    complex*16, dimension(n_basis*(n_basis+1)/2, n_spin), intent(in) :: hamiltonian_w_complex
    integer, intent(in) :: i_band, i_k_point
    real*8, intent(in)  :: i_x,i_y,i_z

!  INPUTS
!   o i_band             -- number of present requested output band
!   o i_k_point          -- number of k-point in present band
!   o i_x, i_y, i_z      -- real(!), relative coordinates of present k-point in units
!                           of the reciprocal lattice vector of the present structure
!  OUTPUT
!    none
!  AUTHOR
!    R.J. Maurer, Yale University 2015, based on output_complex_overlapmatrix.f90 
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE

!  variables needed for HDF5 setup
    integer(HID_T)  :: file_id   ! File identifier for HDF5 usage
    integer(HID_T)  :: plist_id  ! Property list identifier for HDF5 usage
    integer(HID_T)  :: dset_id   ! Dataset identifier
    integer(HID_T)  :: dspace_id ! Dataspace identifier
    integer(HID_T)  :: attr_id   ! Attribute identifier
    integer(HID_T)  :: aspace_id ! Attribute Dataspace identifier
    integer(HID_T)  :: atype_id  ! Attribute Dataspace identifier
    
    integer(HSIZE_T), dimension(5) :: dim_ham ! Dimensions of the hamiltonian matrix
    integer(HSIZE_T), dimension(1) :: adims = (/3/) ! Attribute dimension

    integer         :: drank = 5 ! Dataset rank
    integer         :: arank = 1 ! Attribure rank
    integer(SIZE_T) :: attrlen   ! Length of the attribute string
    character*50    :: file_name ! File name
    character*150, dimension(3) ::  attr_data  ! Attribute data
    character*10    :: num_char
    integer         :: error     ! Error flag

    
!  local variables    
    integer,dimension(:,:),allocatable :: index_b
    complex*16  :: hamiltonian_matrix_entry
    character*150 :: info_str
    !  counters
    integer :: i_row, i_col, i_spin
    integer :: i_entry, i_index
     
!-----------------------------------------------------------------------!
!  begin work
!-----------------------------------------------------------------------!
    if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then

!-----------------------------------------------------------------------!
! Creating the mapping array, only needed for ASCII output
!-----------------------------------------------------------------------!
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

!-----------------------------------------------------------------------!
! Create generic filename for overlap matrix output. Tedious but general.
!-----------------------------------------------------------------------!
      write(file_name, '(A)') 'KS_hamiltonian_matrix.'

! removal of trailing characters.
      write(file_name, '(A)') trim(file_name)

! encode band number
      write(num_char,'(I10)') i_band
      write(file_name, '(A,A,A,A)') trim(file_name),'band_',trim(adjustl(num_char)),'.'

! encode k-point number
      write(num_char,'(I10)') i_k_point
      if (flag_hdf5_unfold) then
        write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.h5'
      else 
        write(file_name, '(A,A,A,A)') trim(file_name),'kpt_',trim(adjustl(num_char)),'.out'
      end if
    
      write (info_str,'(2X,2A)') '| Writing complex hamiltonian matrix: File ', file_name
      call localorb_info ( info_str )
    
!-----------------------------------------------------------------------!
! Start with the writing process
!-----------------------------------------------------------------------!

      if (flag_hdf5_unfold) then

! Define the dimensions of the hamiltonian matrix for HDF5 output
! dim_ham has rank 5, one more than overlap because of spin
        dim_ham(1:5) = (/ n_basis, n_basis, 1, n_spin, 2 /)
! Length of the attribute string
        attrlen = 150

!-----------------------------------------------------------------------!
! Opening H5 data file
!-----------------------------------------------------------------------!
! Initialize HDF5 FORTRAN interface
        call h5open_f(error)
! create the HDF5 file
        call h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
! Create the dataspace.
        call h5screate_simple_f(drank, dim_ham, dspace_id, error)
! Create the dataset with default properties.
        call h5dcreate_f(file_id, 'cplx_hmlt_mtrx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        
!-----------------------------------------------------------------------!
! Create scalar data space for the attribute.
!-----------------------------------------------------------------------!
        call h5screate_simple_f(arank, adims, aspace_id, error)
! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        call h5tset_size_f(atype_id, attrlen, error)
! Create dataset attribute.
        call write_aims_uuid(info_str)
        attr_data(1) = info_str

        write(info_str,'(A,I8)') &
         '# Complex hamiltonian matrix for band number ', i_band
        attr_data(2) = info_str

        write(info_str,'(A,I6,A,3(F11.8,1X))') &
         '# k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', i_x, i_y, i_z
        attr_data(3) = info_str
      
        call h5acreate_f(dset_id, "COM attribute", atype_id, aspace_id, attr_id, error)
! Write the attribute data.  
        call h5awrite_f(attr_id, atype_id, attr_data, adims, error)
! Close the attribute.
        call h5aclose_f(attr_id, error)
!-----------------------------------------------------------------------!

! End access to the dataset and release resources used by it.
        call h5dclose_f(dset_id, error)

        do i_spin = 1, n_spin, 1
            do i_col = 1,n_basis, 1
              do i_row = 1,n_basis, 1
                hamiltonian_matrix_entry = (0.0,0.0)
                i_entry = 0

                if (i_row < i_col) then
                  i_entry = index_b(i_row, i_col)
                  hamiltonian_matrix_entry = hamiltonian_w_complex(i_entry,i_spin)
                else
                  i_entry = index_b(i_col, i_row)
                  hamiltonian_matrix_entry = CONJG(hamiltonian_w_complex(i_entry,i_spin))
                end if
                call out_hmlt_scalapack( &
                     hamiltonian_matrix_entry, &
                     file_id, &
                     plist_id, &
                     i_row, i_col, i_spin, &
                     1, &
                     'cplx_hmlt_mtrx')

              end do ! end loop i_row
            end do ! end loop i_col
        end do ! end loop i_spin
!-----------------------------------------------------------------------!
! Closing H5 data file
!-----------------------------------------------------------------------!
! Terminate access to the data space.
        call h5sclose_f(dspace_id, error)
! Terminate access to the file.
        call h5fclose_f(file_id, error)
! Close FORTRAN interface.!
        call h5close_f(error)

      else
!-----------------------------------------------------------------------!
! File header and initial per-state information first
!-----------------------------------------------------------------------!
        open (15, FILE=file_name)

        write(15,'(A,I8,A,I8,A,3(F12.8,1X))') &
         '# Complex hamiltonian matrix for band number ', i_band, &
         ', k-point number ', i_k_point, ', at relative reciprocal-space coordinates: ', &
         i_x, i_y, i_z
        write(15,'(A)') '#'

!-----------------------------------------------------------------------!
!    find the right entries for every line in the n_basis x n_basis hamiltonian matrix
!-----------------------------------------------------------------------!
        do i_spin = 1, n_spin
            write(15,'(A,I8)') '# spin channel ',i_spin
            do i_row = 1, n_basis
              do i_col = 1, n_basis
                hamiltonian_matrix_entry = (0.0,0.0)
                i_entry = 0
                if (i_row < i_col) then
                  i_entry = index_b(i_row, i_col)
                  hamiltonian_matrix_entry = hamiltonian_w_complex(i_entry,i_spin)
                  write(15,'(3X,F13.8,1X,F13.8)',ADVANCE='NO') &
                            REAL( hamiltonian_matrix_entry ), &
                            AIMAG( hamiltonian_matrix_entry )

                else
                  i_entry = index_b(i_col, i_row)
                  hamiltonian_matrix_entry = hamiltonian_w_complex(i_entry,i_spin)
                  write(15,'(3X,F13.8,1X,F13.8)',ADVANCE='NO') &
                             REAL( hamiltonian_matrix_entry ), &
                             AIMAG(CONJG( hamiltonian_matrix_entry ))
                end if
              end do ! end loop i_col
              write(15,'()') ! end line
            end do ! end loop i_row
        end do ! end loop i_spin
        close(15)
        deallocate(index_b)
      end if ! flag HDF output format
      write (info_str,'(2X,2A)') '| ... finished writing complex hamiltonian matrix: File ', file_name
      call localorb_info ( info_str ) 
    end if ! my_id

  end subroutine output_complex_hamiltonian_matrix
!----------------------------------------------------------------------
!******

!TODO SCALAPACK VERSION IS STILL MISSING

end module hdf5_output
