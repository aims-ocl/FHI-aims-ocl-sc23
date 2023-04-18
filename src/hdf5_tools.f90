!****h* FHI-aims/hdf5_tools
!  NAME
!    hdf5_tools
!  SYNOPSIS


module hdf5_tools
!  PURPOSE
!    Helper subroutines for use when performing parallel output with HDF5
!  AUTHOR
!    The FHI-aims team
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  USES
  use HDF5, only: HID_T, HSIZE_T, SIZE_T ! Data type definitions
  implicit none
!  SOURCE

!  By default, all members below are private
  private

!  Public data type definitions, imported from HDF5 module
  public :: HID_T
  public :: HSIZE_T
  public :: SIZE_T
!  Public subroutines
  public :: open_hdf5
  public :: open_hdf5_dataset
  public :: close_hdf5
  public :: out_k_points_kart
  public :: out_bands_kart
  public :: out_occ_kart
  public :: outmetafile
  public :: open_hdf5_kart
  public :: open_coulelement_kart
  public :: open_coulelement_lvl_v0_1
  public :: out_moment_kart
  public :: out_coulelement_kart
  public :: out_coulelement_lvl_v0_1
  public :: out_KS_vec_scalapack
  public :: out_ovlp_scalapack
  public :: out_hmlt_scalapack
  public :: out_dimensions
  public :: out_basis_atom
  public :: out_KS_eigenvalue

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! General Routines for creating a file and dataset !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine open_hdf5(filename, file_id, plist_id)


  !  PURPOSE
  !   Open hdf5 file
  !
  ! USES
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none


  CHARACTER(*) :: filename  ! File name
  INTEGER(HID_T) :: file_id     ! File identifier
  INTEGER(HID_T) :: plist_id    ! Property list identifier
  ! Open hdf5 file 
  !   
  !  INPUTS
  !    o filename -- file name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

      INTEGER     ::   error  ! Error flag
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, &
                      access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_Independent_F, error)
     

end subroutine open_hdf5

subroutine open_hdf5_dataset(dsetname, file_id,rank, dims)


  !  PURPOSE
  !   Open dataset for Momentummatrix in file with file_id
  !
  ! USES

  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none

  CHARACTER(*) :: dsetname      ! Dataset name
  INTEGER :: rank
  INTEGER(HSIZE_T), DIMENSION(rank)  :: dims
  INTEGER(HID_T) :: file_id                            ! File identifier
  ! Open dataset with rank and dims 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o rank -- rank of dataset
  !    o dims -- dimensions of dataset
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE
      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier


      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL


!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)

end subroutine open_hdf5_dataset

subroutine hdf5_write_attribute(file_id, dsetname, dataset_id, aname,&
                                lines, text)
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none
  CHARACTER(*) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
  INTEGER(HSIZE_T), INTENT(IN) :: lines
  CHARACTER(LEN=200), INTENT(IN) :: text(lines)
  INTEGER(HID_T) :: file_id                            ! File identifier
  CHARACTER(LEN=*) :: aname   ! Attribute name
  ! Write text (dimension(lines)) as atribute to dataset
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o dataset_id -- Dataset identifier
  !    o aname -- Attribute name
  !    o lines -- Number of lines
  !    o text -- Text to write as Attribute to dataset
  !  OUTPUT
  !    o Attribute in Dataset with dsetname and dataset_id in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE
      INTEGER     ::   error  ! Error flag
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      INTEGER(HSIZE_T), DIMENSION(1) :: adims ! Attribute dimension
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attrlen = 200
     adims   = lines
!
!    Initialize FORTRAN interface.
!
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = lines
     CALL h5awrite_f(attr_id, atype_id, text, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
end subroutine hdf5_write_attribute

subroutine close_hdf5(file_id,plist_id  )


  !  PURPOSE
  !   Close hdf5 file
  !
  ! USES
  use HDF5
  implicit none


  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: plist_id  

  ! Close hdf5 file 
  !   
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    
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
  !    Release version, FHI-aims (2008).
  ! SOURCE
     INTEGER     ::   error  ! Error flag
     CALL h5fclose_f(file_id, error)
     CALL h5pclose_f(plist_id, error)
!    Close FORTRAN interface.
!
     CALL h5close_f(error)

end subroutine close_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Routines for writting k_points_list and E_bands !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine out_k_points(file_id,plist_id)

  !  PURPOSE
  !   Write k-points list to file
  !
  ! USES

  use dimensions, only: n_k_points
  use HDF5
  use pbc_lists, only: k_point_list
  implicit none

  !  ARGUMENTS
  ! Write k-Points list to file
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o file dipelement_k_point.dat
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     INTEGER(HID_T) :: plist_id   ! Property list identifier
     CHARACTER(LEN=8) :: dsetname = "k_points"      ! Dataset name
     INTEGER(HID_T) :: file_id        ! File identifier
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k
     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(2) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 2 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
     allocate(dset1_data(4,n_k_points))
     dims1 = (/4,n_k_points/)
     do i_k = 1, n_k_points
          dset1_data(1,i_k)= i_k
          dset1_data(2,i_k)= k_point_list(i_k,1)
          dset1_data(3,i_k)= k_point_list(i_k,2)
          dset1_data(4,i_k)= k_point_list(i_k,3)
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = 4
     data_dims(2) = n_k_points
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_k_points

subroutine out_k_points_kart(file_id,plist_id)

  !  PURPOSE
  !   Write k-points list to file
  !
  ! USES
  use HDF5
  use pbc_lists, only: k_point_list
  use runtime_choices, only: n_k_points_xyz
  implicit none

  !  ARGUMENTS
  ! Write k-Points list to file
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o file dipelement_k_point.dat
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     INTEGER(HID_T) :: plist_id   ! Property list identifier
     CHARACTER(LEN=8) :: dsetname = "k_points"      ! Dataset name
     INTEGER(HID_T) :: file_id        ! File identifier
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k, i_x, i_y, i_z
     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:,:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(4) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 4 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(4) :: data_dims

     INTEGER(HID_T) :: attr_id       ! Attribute identifier
     INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
     INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/2/) ! Attribute dimension
     INTEGER     ::   arank = 1                      ! Attribure rank
     INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
     CHARACTER(LEN=120), DIMENSION(2) ::  attr_data  ! Attribute data
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
     CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name

     attr_data(1) = " k_x x k_y x k_z x (i_k, k_x, k_y, k_z)"
     attr_data(2) = "k_grid in relative units of the reziprocal lattice"//&
                    " vectors)"
     attrlen = 120
     allocate(dset1_data(4,n_k_points_xyz(1),n_k_points_xyz(2),&
              n_k_points_xyz(3)))
     dims1 = (/4,n_k_points_xyz(1),n_k_points_xyz(2),n_k_points_xyz(3)/)
     i_k=0
     do i_x = 1, n_k_points_xyz(1)
        do i_y = 1, n_k_points_xyz(2)
           do i_z = 1, n_k_points_xyz(3) 
              i_k=i_k+1  
              dset1_data(1,i_x,i_y,i_z)=i_k-1
              dset1_data(2,i_x,i_y,i_z)=k_point_list(i_k,1)
              dset1_data(3,i_x,i_y,i_z)=k_point_list(i_k,2)
              dset1_data(4,i_x,i_y,i_z)=k_point_list(i_k,3)
           end do
        end do
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = 4
     data_dims(2) = n_k_points_xyz(1)
     data_dims(3) = n_k_points_xyz(2)
     data_dims(4) = n_k_points_xyz(3)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
			      attr_id, error)
     data_dims_attr(1) = 2
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_k_points_kart

subroutine out_bands(KS_Eigenvalue,file_id,plist_id, n_state_min_in, &
                     n_state_max_in)

  !  PURPOSE
  !   Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !
  ! USES
  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  use pbc_lists, only: k_point_list
  implicit none

  !  ARGUMENTS
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::  KS_eigenvalue
  ! Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !  INPUTS
  !    o KS_Eigenvalue -- KS-Eigenvalues
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     CHARACTER(LEN=7) :: dsetname = "E_bands"      ! Dataset name
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k

     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(2) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 2 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

     allocate(dset1_data(4+(n_state_max_in-n_state_min_in+1)+1,n_k_points))
     dims1 = (/4+(n_state_max_in-n_state_min_in+1)+1,n_k_points/)
     do i_k = 1, n_k_points
          dset1_data(1,i_k)= i_k
          dset1_data(2,i_k)= k_point_list(i_k,1)
          dset1_data(3,i_k)= k_point_list(i_k,2)
          dset1_data(4,i_k)= k_point_list(i_k,3)
          dset1_data(5:(n_state_max_in-n_state_min_in+1)+1+5,i_k)=&
                        KS_eigenvalue(n_state_min_in:n_state_max_in, 1, i_k)
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = 4+(n_state_max_in-n_state_min_in+1)+1
     data_dims(2) = n_k_points
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error,xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_bands

subroutine out_occ(occ_numbers,file_id,plist_id, n_state_min_in, &
                     n_state_max_in)

  !  PURPOSE
  !   Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !
  ! USES

  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  use pbc_lists, only: k_point_list
  implicit none

  !  ARGUMENTS
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  real*8, dimension(n_states, n_spin,n_k_points) :: occ_numbers
  ! Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !  INPUTS
  !    o occ_numbers -- occ_numbers
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     CHARACTER(LEN=18) :: dsetname = "occupation_numbers"      ! Dataset name
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k

     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(2) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 2 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

     allocate(dset1_data(4+(n_state_max_in-n_state_min_in+1)+1,n_k_points))
     dims1 = (/4+(n_state_max_in-n_state_min_in+1)+1,n_k_points/)
     do i_k = 1, n_k_points
          dset1_data(1,i_k)= i_k
          dset1_data(2,i_k)= k_point_list(i_k,1)
          dset1_data(3,i_k)= k_point_list(i_k,2)
          dset1_data(4,i_k)= k_point_list(i_k,3)
          dset1_data(5:(n_state_max_in-n_state_min_in+1)+1+5,i_k)=&
                        occ_numbers(n_state_min_in:n_state_max_in, 1, i_k)
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = 4+(n_state_max_in-n_state_min_in+1)+1
     data_dims(2) = n_k_points
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error,xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_occ

subroutine out_bands_kart(KS_Eigenvalue,file_id,plist_id,n_state_min_in,&
                          n_state_max_in)

  !  PURPOSE
  !   Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !
  ! USES

  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  use runtime_choices, only: n_k_points_xyz
  implicit none

  !  ARGUMENTS
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::  KS_eigenvalue
  ! Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !  INPUTS
  !    o KS_Eigenvalue -- KS-Eigenvalues
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     CHARACTER(LEN=7) :: dsetname = "E_bands"      ! Dataset name
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k, i_x, i_y, i_z, i_spin

     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HID_T) :: attr_id       ! Attribute identifier
     INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
     INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/5/) ! Attribute dimension
     INTEGER     ::   arank = 1                      ! Attribure rank
     INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
     CHARACTER(LEN=120), DIMENSION(5) ::  attr_data  ! Attribute data
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
     CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name

     attr_data(1) = "spin addition"
     attr_data(2) = "k_x x k_y x k_z x Energy Eigenvalues in [Ha]"
     attr_data(3) = "see 'k_points' dataset for k_grid (in relative units of"//&
                    " the reziprocal lattice vectors)"
     attr_data(4) = "see 'CBM' dataset for CBM bandindex"
     attr_data(5) = "see 'VBM' dataset for VBM bandindex"
     attrlen = 120

     allocate(dset1_data(n_spin, (n_state_max_in-n_state_min_in+1),n_k_points_xyz(1),&
                          n_k_points_xyz(2),n_k_points_xyz(3)))
     dims1 = (/n_spin, (n_state_max_in-n_state_min_in+1),n_k_points_xyz(1),&
                n_k_points_xyz(2),n_k_points_xyz(3)/)
     i_k=0
     do i_x = 1, n_k_points_xyz(1)
        do i_y = 1, n_k_points_xyz(2)
           do i_z = 1, n_k_points_xyz(3) 
              i_k=i_k+1  
             do i_spin = 1, n_spin
              dset1_data(i_spin,:,i_x,i_y,i_z)=KS_eigenvalue(n_state_min_in:&
                                                    n_state_max_in, i_spin, i_k)
             end do
           end do
        end do
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = n_spin
     data_dims(2) = (n_state_max_in-n_state_min_in+1)
     data_dims(3) = n_k_points_xyz(1)
     data_dims(4) = n_k_points_xyz(2)
     data_dims(5) = n_k_points_xyz(3)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims,&
                     error, xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
			      attr_id, error)
     data_dims_attr(1) = 5
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_bands_kart

subroutine out_occ_kart(occ_numbers,file_id,plist_id,n_state_min_in,&
                          n_state_max_in)

  !  PURPOSE
  !   Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !
  ! USES

  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  use runtime_choices, only: n_k_points_xyz
  implicit none

  !  ARGUMENTS
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  ! Write band structure to file for all k-points (n_state_min to 
  !   n_state_max)
  !  INPUTS
  !    o KS_Eigenvalue -- KS-Eigenvalues
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     CHARACTER(LEN=18) :: dsetname = "occupation_numbers"      ! Dataset name
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k, i_x, i_y, i_z, i_spin

     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HID_T) :: attr_id       ! Attribute identifier
     INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
     INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/5/) ! Attribute dimension
     INTEGER     ::   arank = 1                      ! Attribure rank
     INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
     CHARACTER(LEN=120), DIMENSION(5) ::  attr_data  ! Attribute data
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
     CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name

     attr_data(1) = "spin addition"
     attr_data(1) = "k_x x k_y x k_z x Energy Eigenvalues in [Ha]"
     attr_data(2) = "see 'k_points' dataset for k_grid (in relative units of"//&
                    " the reziprocal lattice vectors)"
     attr_data(3) = "see 'CBM' dataset for CBM bandindex"
     attr_data(4) = "see 'VBM' dataset for VBM bandindex"
     attrlen = 120

     allocate(dset1_data(n_spin, (n_state_max_in-n_state_min_in+1),n_k_points_xyz(1),&
                          n_k_points_xyz(2),n_k_points_xyz(3)))
     dims1 = (/n_spin, (n_state_max_in-n_state_min_in+1),n_k_points_xyz(1),&
                n_k_points_xyz(2),n_k_points_xyz(3)/)
     i_k=0
     do i_x = 1, n_k_points_xyz(1)
        do i_y = 1, n_k_points_xyz(2)
           do i_z = 1, n_k_points_xyz(3) 
              i_k=i_k+1  
              do i_spin = 1, n_spin
                dset1_data(i_spin,:,i_x,i_y,i_z)=occ_numbers(n_state_min_in:&
                                                      n_state_max_in, i_spin, i_k)
              end do
           end do
        end do
     end do
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = n_spin
     data_dims(2) = (n_state_max_in-n_state_min_in+1)
     data_dims(3) = n_k_points_xyz(1)
     data_dims(4) = n_k_points_xyz(2)
     data_dims(5) = n_k_points_xyz(3)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims,&
                     error, xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
			      attr_id, error)
     data_dims_attr(1) = 5
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_occ_kart

subroutine outmetafile(KS_eigenvalue,occ_numbers,chemical_potential,file_id,&
                       plist_id,n_state_min_in, n_state_max_in)

  !  PURPOSE
  !   Write metadata to hdf5 file
  !
  ! USES
  use constants, only: bohr
  use dimensions, only: n_states, n_spin, n_k_points, n_periodic
  use geometry, only: recip_lattice_vector
  use HDF5
  implicit none

  integer :: n_state_min_in
  integer :: n_state_max_in
  real*8, dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue 
  real*8, dimension(n_states, n_spin,n_k_points) :: occ_numbers
   real*8 :: chemical_potential
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  INTEGER(HID_T) :: file_id        ! File identifier
  ! Close hdf5 file 
  !   
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o chemical_potential -- chemical potential
  !    o KS_eigenvalue -- KS eigenvalues
  !    o occ_numbers -- occupation numbers
  !  OUTPUT
  !    
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


      character*150 :: info_str
      real*8, dimension(n_spin) :: n_elec_per_channel

      real*8 :: homo_level, lumo_level
      real*8 :: homo_occ, lumo_occ
      integer :: i_kpt_homo, i_kpt_lumo
      real*8 :: i_state_homo, i_state_lumo
      integer :: i_spin_homo, i_spin_lumo
      real*8 :: midpoint

      integer i_state
      integer i_spin
      integer i_k_point
 
     CHARACTER(LEN=3) :: dsetname1 = "CBM"      ! Dataset name
     CHARACTER(LEN=3) :: dsetname2 = "VBM"      ! Dataset name
     CHARACTER(LEN=12) :: dsetname3 = "Fermi_energy"      ! Dataset name
     CHARACTER(LEN=13) :: dsetname4 = "Energy_window"      ! Dataset name
     CHARACTER(LEN=28) :: dsetname6 = "reziprocal_lattice_vectors"! Dataset name

     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k

     INTEGER     ::   error ! Error flag
     REAL*8, DIMENSION(3) :: dset1_data  ! Arrays to hold data
     REAL*8, DIMENSION(3) :: dset2_data
     REAL*8, DIMENSION(1) :: dset3_data
     REAL*8, DIMENSION(1,2) :: dset4_data
     REAL*8, DIMENSION(3,3) :: dset6_data
     INTEGER(HSIZE_T), DIMENSION(1) :: dims1 ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(1) :: dims2 ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(1) :: dims3 ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(2) :: dims4 ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(2) :: dims6 ! Dataset dimensions
     INTEGER     ::   rank1 = 1 ! Datasets rank
     INTEGER     ::   rank2 = 2 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(1) :: data1_dims
     INTEGER(HSIZE_T), DIMENSION(1) :: data2_dims
     INTEGER(HSIZE_T), DIMENSION(1) :: data3_dims
     INTEGER(HSIZE_T), DIMENSION(2) :: data4_dims
     INTEGER(HSIZE_T), DIMENSION(2) :: data6_dims

     INTEGER(HID_T) :: attr_id       ! Attribute identifier
     INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
     INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/3/) ! Attribute dimension
     INTEGER     ::   arank = 1                      ! Attribure rank
     INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
     CHARACTER(LEN=120), DIMENSION(3) ::  attr_data  ! Attribute data
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
     CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name


      homo_level = -10000.
      lumo_level = 10000.

      ! Define the correct "half occupation" (with or without spin)
      if (n_spin.eq.1) then
        midpoint = 1.d0
      else
        midpoint = 0.5d0
      end if
      i_kpt_homo = 0
      i_kpt_lumo = 0
      i_state_homo = 0
      i_state_lumo = 0
      do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1

            if (occ_numbers(i_state,i_spin,i_k_point).ge.midpoint) then
              ! check if homo
              ! "HOMO" also includes Fermi level ("ge" above)
              if (KS_eigenvalue(i_state,i_spin,i_k_point).gt.homo_level)&
                    then
                homo_level = KS_eigenvalue(i_state,i_spin,i_k_point)
                homo_occ = occ_numbers(i_state,i_spin,i_k_point)
                i_kpt_homo = i_k_point
                i_spin_homo = i_spin
                i_state_homo = i_state
              end if
            else
            ! check if lumo
              if (KS_eigenvalue(i_state,i_spin,i_k_point).lt.lumo_level)&
                     then
                lumo_level = KS_eigenvalue(i_state,i_spin,i_k_point)
                lumo_occ = occ_numbers(i_state,i_spin,i_k_point)
                i_kpt_lumo = i_k_point
                i_spin_lumo = i_spin
                i_state_lumo = i_state
              end if
            end if

          enddo
        enddo
      enddo

     attr_data(1) = "0: CBM energy in [Ha] "
     attr_data(2) = "1: CBM @ k_point number (see 'k_points' dataset for"//&
                                            " relative coordinates)"
     attr_data(3) = "2: Bandindex of CBM with respect to 'E_bands'"
     attrlen = 120

     dims1 = (/3/)
     dset1_data(1)= lumo_level
     dset1_data(2)= i_kpt_lumo
     dset1_data(3)= i_state_lumo-n_state_min_in+1
     CALL h5screate_simple_f(rank1, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data1_dims(1) = 3
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data1_dims,&
                     error,xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 3
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)

     attr_data(1) = "0: VBM energy in [Ha] "
     attr_data(2) = "1: VBM @ k_point number (see 'k_points' dataset for"//&
                                            " relative coordinates)"
     attr_data(3) = "2: Bandindex of VBM with respect to 'E_bands'"
     attrlen = 120
     dims2 = (/3/)
     dset2_data(1)= homo_level
     dset2_data(2)= i_kpt_homo
     dset2_data(3)= i_state_homo-n_state_min_in+1
     CALL h5screate_simple_f(rank1, dims2, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data2_dims(1) = 3
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset2_data, data2_dims,&
                     error,xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 3
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)

     attr_data(1) = "Fermi energy in [Ha]"
     attr_data(2) = " "
     attr_data(3) = " "
     attrlen = 120
     dims3 = (/1/)
     dset3_data(1)= chemical_potential
     CALL h5screate_simple_f(rank1, dims3, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data3_dims(1) = 1
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset3_data, data3_dims,&
                     error,xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 3
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)

     attr_data(1) = "Preselected Energy_window"
     attr_data(2) = "0: Lower bound i state"
     attr_data(3) = "1: Upper bound j_state"
     attrlen = 120
     dims4 = (/1,2/)
     dset4_data(1,1)= n_state_min_in
     dset4_data(1,2)= n_state_max_in
     CALL h5screate_simple_f(rank2, dims4, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data4_dims(1) = 1
     data4_dims(1) = 2
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset4_data, data4_dims,&
                     error, xfer_prp = plist_id)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 3
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)

     if (n_periodic .gt. 0) then
	attr_data(1) = "0: k_x in inverse Angstrom"
	attr_data(2) = "1: k_y in inverse Angstrom"
	attr_data(3) = "2: k_z in inverse Angstrom"
	attrlen = 120
	dims6 = (/3,3/)
	dset6_data(:,1)= recip_lattice_vector(:,1)/bohr
	dset6_data(:,2)= recip_lattice_vector(:,2)/bohr
	dset6_data(:,3)= recip_lattice_vector(:,3)/bohr
	CALL h5screate_simple_f(rank2, dims6, dataspace_id, error)
	CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dataspace_id, &
				  dataset_id, error)
	data6_dims(1) = 3
	data6_dims(2) = 3
	CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset6_data, data6_dims,&
			error, xfer_prp = plist_id)
	CALL h5screate_simple_f(arank, adims, aspace_id, error)
	CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
	CALL h5tset_size_f(atype_id, attrlen, error)
	CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
			  attr_id, error)
	data_dims_attr(1) = 3
	CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
	CALL h5aclose_f(attr_id, error)
	CALL h5sclose_f(aspace_id, error)
	CALL h5sclose_f(dataspace_id, error)
	CALL h5dclose_f(dataset_id, error)
      endif

end subroutine outmetafile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Specialized routine for writing Momentummatrix !!!!!!!!!!!!!!!!!!
!!!!!!!!!!! and Dipolematrix to datasets                   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine open_hdf5_state(dsetname, file_id, plist_id,n_state_min_in,&
                           n_state_max_in)


  !  PURPOSE
  !   Open dataset for Momentummatrix in file with file_id
  !
  ! USES

  use dimensions, only: n_spin, n_k_points
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none

  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER, INTENT(IN) :: n_state_min_in
  INTEGER, INTENT(IN) :: n_state_max_in
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  ! Open dataset for Momentummatrix in file with file_id 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE
      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier

      INTEGER(HSIZE_T), DIMENSION(3) :: dims ! Dataset dimensions
      INTEGER     ::   rank = 3 ! Datasets rank
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/6/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      CHARACTER(LEN=150), DIMENSION(6) ::  attr_data  ! Attribute data
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attr_data(1) = "Dipolematrix/Momentummatrix"
     attr_data(2) = "k_x"
     attr_data(3) = "k_y"
     attr_data(4) = "k_z"
     attr_data(5) = "E_i -> E_j (upper Triangle of matrix: num=0;"//&
                            "do i_spin=1, n_spin; do j_spin=i_spin, n_spin;"//&
                            "do i=0, N-1; do j=i, N-1;p(num)=X(i,j);"//&
                            "num=num+1; enddo;enddo; enddo;enddo)"
     attr_data(6) = "(Re{Momentummatrix_x},Im{Momentummatrix_x},"//&
                    "Re{Momentummatrix_y},Im{Momentummatrix_y},"//&
                    "Re{Momentummatrix_z},Im{Momentummatrix_z})"
     attrlen = 200

!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     dims = (/6,(((n_state_max_in-n_state_min_in+1)+1)*&
                  (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)&
                  ,n_k_points/)
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
			      attr_id, error)
     data_dims_attr(1) = 6
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)

end subroutine open_hdf5_state

subroutine open_hdf5_kart(dsetname, file_id, plist_id,n_state_min_in,&
                          n_state_max_in)


  !  PURPOSE
  !   Open dataset for Momentummatrix in file with file_id
  !
  ! USES

  use dimensions, only: n_spin
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  use runtime_choices, only: n_k_points_xyz
  implicit none

  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER, INTENT(IN) :: n_state_min_in
  INTEGER, INTENT(IN) :: n_state_max_in
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  ! Open dataset for Momentummatrix in file with file_id 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_in -- Minimum state
  !    o n_state_max_in -- Maximum state
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE
      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      INTEGER(HSIZE_T), DIMENSION(5) :: dims ! Dataset dimensions
      INTEGER     ::   rank = 5 ! Datasets rank
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/6/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      CHARACTER(LEN=150), DIMENSION(6) ::  attr_data  ! Attribute data
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attr_data(1) = "Dipolematrix/Momentummatrix"
     attr_data(2) = "k_x"
     attr_data(3) = "k_y"
     attr_data(4) = "k_z"
     attr_data(5) = "E_i -> E_j (upper Triangle of matrix: num=0;"//&
                            "do i_spin=1, n_spin; do j_spin=i_spin, n_spin;"//&
                            "do i=0, N-1; do j=i, N-1;p(num)=X(i,j);"//&
                            "num=num+1; enddo;enddo; enddo;enddo)"
     attr_data(6) = "(Re{Momentummatrix_x},Im{Momentummatrix_x},"//&
                     "Re{Momentummatrix_y},Im{Momentummatrix_y},"//&
                     "Re{Momentummatrix_z},Im{Momentummatrix_z})"
     attrlen = 200

!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     dims = (/6,(((n_state_max_in-n_state_min_in+1)+1)*&
                  (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),&
                   n_k_points_xyz(1),n_k_points_xyz(2),n_k_points_xyz(3)/)
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 6
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)
end subroutine open_hdf5_kart

subroutine open_coulelement_kart(dsetname, file_id, plist_id, &
                                 n_state_min_coulmat, n_state_max_coulmat, &
                                 n_q_point_in)


  !  PURPOSE
  !   Open dataset for Fourier component Coulomb matrix in file with file_id
  !
  ! USES

  use dimensions, only: n_k_points
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none

      CHARACTER(*) :: dsetname      ! Dataset name
      INTEGER(HID_T) :: file_id                            ! File identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      INTEGER, INTENT(IN) :: n_state_min_coulmat
      INTEGER, INTENT(IN) :: n_state_max_coulmat
      INTEGER, INTENT(IN) :: n_q_point_in

  ! Open dataset for Momentummatrix in file with file_id 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_coulmat -- Minimum state
  !    o n_state_max_coulmat -- Maximum state
  !    o n_q_point_in -- Number of q-points
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      INTEGER(HSIZE_T), DIMENSION(5) :: dims ! Dataset dimensions
      INTEGER     ::   rank = 5 ! Datasets rank
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/6/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      CHARACTER(LEN=150), DIMENSION(6) ::  attr_data  ! Attribute data
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attr_data(1) = "Coulombmatrix"
     attr_data(2) = "n_q"
     attr_data(3) = "n_k"
     attr_data(4) = "n_k_strich"
     attr_data(5) = "do n_state=n_state_min,n_state_max;do m_state=n_state, "//&
                    "n_state_max enddo;enddo"
     attr_data(6) = "(Re{Coulombmatrix},Im{Coulombmatrix})"
     attrlen = 150

!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     dims = (/2,((n_state_max_coulmat-n_state_min_coulmat+1)+1)*&
            (n_state_max_coulmat-n_state_min_coulmat+1)/2,n_q_point_in,&
             n_k_points,n_k_points/)
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 6
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)

end subroutine open_coulelement_kart

subroutine open_coulelement_lvl_v0(dsetname, file_id, plist_id, n_q_point_in&
                                   )


  !  PURPOSE
  !   Open dataset for Fourier component Coulomb matrix in file with file_id
  !
  ! USES
  use dimensions, only: n_k_points
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none

      CHARACTER(LEN=14) :: dsetname      ! Dataset name
      INTEGER(HID_T) :: file_id                            ! File identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      INTEGER, INTENT(IN) :: n_q_point_in

  ! Open dataset for Coulombmatrix (from product basis) in file with file_id 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_q_point_in -- Number of q-points
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      INTEGER(HSIZE_T), DIMENSION(4) :: dims ! Dataset dimensions
      INTEGER     ::   rank = 4 ! Datasets rank
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/5/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      CHARACTER(LEN=150), DIMENSION(6) ::  attr_data  ! Attribute data
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attr_data(1) = "Coulombmatrix"
     attr_data(2) = "n_q"
     attr_data(3) = "n_k"
     attr_data(4) = "n_k_strich"
     attr_data(5) = "(Re{Coulombmatrix},Im{Coulombmatrix})"
     attrlen = 150

!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     dims = (/2,n_q_point_in,n_k_points,n_k_points/)
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 5
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)

end subroutine open_coulelement_lvl_v0

subroutine open_coulelement_lvl_v0_1(dsetname, file_id, plist_id, n_q_point_in,&
                                   size_KS)


  !  PURPOSE
  !   Open dataset for Fourier component Coulomb matrix in file with file_id
  !
  ! USES

  use dimensions, only: n_k_points
  use HDF5
  use mpi_tasks, only: mpi_comm_world, mpi_info_null
  implicit none

      CHARACTER(LEN=14) :: dsetname      ! Dataset name
      INTEGER(HID_T) :: file_id                            ! File identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      INTEGER, INTENT(IN) :: n_q_point_in
      INTEGER, INTENT(IN) :: size_KS

  ! Open dataset for Coulombmatrix (from product basis) in file with file_id 
  !   
  !  INPUTS
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_q_point_in -- Number of q-points
  !  OUTPUT
  !    o Dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

      INTEGER     ::   error  ! Error flag
      INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      INTEGER(HSIZE_T), DIMENSION(5) :: dims ! Dataset dimensions
      INTEGER     ::   rank = 5 ! Datasets rank
      INTEGER :: mpierror       ! MPI error flag
      INTEGER :: comm, info
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/6/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      CHARACTER(LEN=150), DIMENSION(9) ::  attr_data  ! Attribute data
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims_attr
      CHARACTER(LEN=8), PARAMETER :: aname = "contents"   ! Attribute name
     !
     ! Initialize FORTRAN interface.
     !
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL

     attr_data(1) = "Coulombmatrix"
     attr_data(2) = "n_q"
     attr_data(3) = "n_k"
     attr_data(4) = "n_k_strich"
     attr_data(5) = "n_state_min:n_state_max**4"
     attr_data(6) = "(Re{Coulombmatrix},Im{Coulombmatrix})"
     attrlen = 150

!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     !
     dims = (/2,size_KS*size_KS*size_KS*size_KS,n_q_point_in,n_k_points,&
              n_k_points/)
     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     CALL h5screate_simple_f(arank, adims, aspace_id, error)
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
     CALL h5tset_size_f(atype_id, attrlen, error)
     CALL h5acreate_f(dataset_id, aname, atype_id, aspace_id, &
		      attr_id, error)
     data_dims_attr(1) = 9
     CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims_attr, error)
     CALL h5aclose_f(attr_id, error)
     CALL h5sclose_f(aspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     CALL h5sclose_f(dataspace_id, error)

end subroutine open_coulelement_lvl_v0_1

subroutine out_moment_state(dipelement,KS_Eigenvalue, k_point,dsetname,&
                                file_id,plist_id,i_coord,n_state_min_in,&
                                n_state_max_in)


  !  PURPOSE
  !  Write Momentummatrix elements to dataset in file for one k-point (parallel)
  !
  ! USES

  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  implicit none

  !  ARGUMENTS
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer, INTENT(IN) :: i_coord
  complex*16, dimension((((n_state_max_in-n_state_min_in+1)+1)*&
                (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),&
                INTENT(IN) ::  dipelement
  integer, INTENT(IN) :: k_point
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::  KS_eigenvalue
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  CHARACTER(LEN=14) :: dsetname     ! Dataset name
  ! Write Momentummatrix-elements to file
  !  INPUTS
  !    o dipelement_one -- (n_states+1)*n_states/2 Momentummatrix (x)
  !    o dipelement_two -- (n_states+1)*n_states/2 Momentummatrix (y)
  !    o dipelement_three -- (n_states+1)*n_states/2 Momentummatrix (z)
  !    o KS_eigenvalue -- KS-Eigenvalues
  !    o k_point -- k-point to be written
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o dsetname -- Dataset name
  !  OUTPUT
  !    o Dataset dsetname in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: n_state, m_state, num, ierror

!HDF5 infrastructure

     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::  i, j

     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank
     REAL*8, ALLOCATABLE  :: dset1_data(:,:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(3) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 3 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(3) :: offset

     dims1 = (/2,(((n_state_max_in-n_state_min_in+1)+1)*(n_state_max_in-&
               n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),1/)
     allocate (dset1_data(2,(((n_state_max_in-n_state_min_in+1)+1)*&
            (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),1))
     num = 0
     do n_state = n_state_min_in,n_state_max_in
       do m_state = n_state, n_state_max_in
          num = num + 1
          dset1_data(1,num,1)= Real(dipelement(num))
          dset1_data(2,num,1)= aImag(dipelement(num))
      end do
     end do

     CALL h5dopen_f(file_id, dsetname, dataset_id, error)

     data_dims(1) = 2
     data_dims(2) = (((n_state_max_in-n_state_min_in+1)+1)*&
                    (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)
     data_dims(3) = 1
     offset(1) = 2*(i_coord-1)
     offset(2) = 0
     offset(3) = k_point-1
     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims,&
                      error,file_space_id = dataspace_id, &
                      mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)

     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_moment_state

subroutine out_moment_kart(dipelement, KS_Eigenvalue, k_point,dsetname,&
                           file_id,plist_id,coord,n_state_min_in,n_state_max_in)

  !  PURPOSE
  !  Write Momentummatrix elements to dataset in file for one k-point (parallel)
  !
  ! USES

  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  use pbc_lists, only: k_point_list
  use runtime_choices, only: n_k_points_xyz
  implicit none

  !  ARGUMENTS
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  complex*16, dimension((((n_state_max_in-n_state_min_in+1)+1)*&
              (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)), &
                                                    INTENT(IN) ::     dipelement
  integer, INTENT(IN) :: k_point
  integer, INTENT(IN) :: coord
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::  KS_eigenvalue
     
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  CHARACTER(LEN=14) :: dsetname      ! Dataset name

 ! Write Momentummatrix-elements to file
  !  INPUTS
  !    o dipelement_one -- (n_states+1)*n_states/2 Momentummatrix (x)
  !    o dipelement_two -- (n_states+1)*n_states/2 Momentummatrix (y)
  !    o dipelement_three -- (n_states+1)*n_states/2 Momentummatrix (z)
  !    o KS_eigenvalue -- KS-Eigenvalues
  !    o k_point -- k-point to be written
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o dsetname -- Dataset name
  !  OUTPUT
  !    o Dataset dsetname in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     integer:: n_state, m_state, i_spin, j_spin, num, ierror


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::  i, j

     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank
     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data



     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(5) :: offset

     dims1 = (/2,((((n_state_max_in-n_state_min_in+1)+1)*&
           (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),1,1,1/)
     allocate(dset1_data(2,((((n_state_max_in-n_state_min_in+1)+1)*&
            (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),1,1,1))
     num = 0
     do i_spin = 1, n_spin
       do j_spin = i_spin , n_spin
	do n_state = n_state_min_in,n_state_max_in
	  do m_state = n_state, n_state_max_in
	      num = num + 1
	      dset1_data(1,num,1,1,1)= Real(dipelement(num))
	      dset1_data(2,num,1,1,1)= aImag(dipelement(num))
	  end do
	end do
       enddo
     enddo

     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 2
     data_dims(2) = (((n_state_max_in-n_state_min_in+1)+1)*&
                    (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)
     data_dims(3) = 1
     data_dims(4) = 1
     data_dims(5) = 1
     offset(1) = (coord-1)*2
     offset(2) = 0
     offset(3) = anint((k_point_list(k_point,1)-k_point_list(1,1))*&
                                                          n_k_points_xyz(1))
     offset(4) = anint((k_point_list(k_point,2)-k_point_list(1,2))*&
                                                          n_k_points_xyz(2))
     offset(5) = anint((k_point_list(k_point,3)-k_point_list(1,3))*&
                                                          n_k_points_xyz(3))
     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_moment_kart

subroutine out_coulelement_kart( coul_val, q_point, k_point, kstrich_point,&
                                dsetname,file_id,plist_id, n_state_min_coulmat,&
                                n_state_max_coulmat)

  !  PURPOSE
  !   Write Fourier component of Coulombelement to file for one q-point, 
  !   k-point, kstrich_point
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  complex*16, INTENT(IN) :: coul_val((&
                                (n_state_max_coulmat-n_state_min_coulmat+1)+1)&
                               *(n_state_max_coulmat-n_state_min_coulmat+1)/2)
  integer, INTENT(IN) :: q_point
  integer, INTENT(IN) :: k_point
  integer, INTENT(IN) :: kstrich_point
  CHARACTER(*) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  INTEGER, INTENT(IN) :: n_state_min_coulmat
  INTEGER, INTENT(IN) :: n_state_max_coulmat
  ! Write dipelement to file
  !  INPUTS
  !    o coul_val -- coulomb element at q, k, k'
  !    o q_point -- q-point
  !    o k_point -- k-point
  !    o kstrich_point -- kstrich-point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !    o n_state_min_coulmat -- Minimum state
  !    o n_state_max_coulmat -- Maximum state
  !  OUTPUT
  !    o Data in dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(5) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     INTEGER :: n_state, m_state, num
     num=0
     dims1 = (/2,((n_state_max_coulmat-n_state_min_coulmat+1)+1)*&
                  (n_state_max_coulmat-n_state_min_coulmat+1)/2,1,1,1/)
     allocate(dset1_data(2,((n_state_max_coulmat-n_state_min_coulmat+1)+1)*&
                           (n_state_max_coulmat-n_state_min_coulmat+1)/2,1,1,1))
     do n_state = n_state_min_coulmat, n_state_max_coulmat
       do m_state = n_state, n_state_max_coulmat
	 num = num + 1
         dset1_data(1,num,1,1,1)=Real(coul_val(num))
         dset1_data(2,num,1,1,1)=Aimag(coul_val(num))
       enddo
     enddo  
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 2
     data_dims(2) = ((n_state_max_coulmat-n_state_min_coulmat+1)+1)*&
                     (n_state_max_coulmat-n_state_min_coulmat+1)/2
     data_dims(3) = 1
     data_dims(4) = 1
     data_dims(5) = 1
     offset(1) = 0
     offset(2) = 0
     offset(3) = q_point-1
     offset(4) = k_point-1
     offset(5) = kstrich_point-1

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_coulelement_kart

subroutine out_coulelement_lvl_v0( coul_val, q_point, k_point, kstrich_point,&
                                dsetname,file_id,plist_id)

  !  PURPOSE
  !   Coulombelement (from lvl method) to file for one q-point, 
  !   k-point, kstrich_point
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  complex*16, INTENT(IN) :: coul_val
  integer, INTENT(IN) :: q_point
  integer, INTENT(IN) :: k_point
  integer, INTENT(IN) :: kstrich_point
  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  ! Write dipelement to file
  !  INPUTS
  !    o coul_val -- coulomb element at q, k, k'
  !    o q_point -- q-point
  !    o k_point -- k-point
  !    o kstrich_point -- kstrich-point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Data in dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(4) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 4 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(4) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(4) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:)  ! Arrays to hold data

     dims1 = (/2,1,1,1/)
     allocate(dset1_data(2,1,1,1))
     dset1_data(1,1,1,1)=Real(coul_val)
     dset1_data(2,1,1,1)=Aimag(coul_val)
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 2
     data_dims(2) = 1
     data_dims(3) = 1
     data_dims(4) = 1
     offset(1) = 0
     offset(2) = q_point-1
     offset(3) = k_point-1
     offset(4) = kstrich_point-1

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_coulelement_lvl_v0

subroutine out_coulelement_lvl_v0_1( coul_val, q_point, k_point, kstrich_point,&
                                dsetname,file_id,plist_id,size_KS)

  !  PURPOSE
  !   Coulombelement (from lvl method) to file for one q-point, 
  !   k-point, kstrich_point
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  integer, INTENT(IN) :: size_KS
  complex*16, INTENT(IN) :: coul_val(size_KS,size_KS,size_KS,size_KS)
  integer, INTENT(IN) :: q_point
  integer, INTENT(IN) :: k_point
  integer, INTENT(IN) :: kstrich_point
  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  ! Write dipelement to file
  !  INPUTS
  !    o coul_val -- coulomb element at q, k, k'
  !    o q_point -- q-point
  !    o k_point -- k-point
  !    o kstrich_point -- kstrich-point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Data in dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(5) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     INTEGER     ::   i_state, j_state, k_state, l_state, num
     dims1 = (/2,size_KS*size_KS*size_KS*size_KS,1,1,1/)
     allocate(dset1_data(2,size_KS*size_KS*size_KS*size_KS,1,1,1))
     num=0
     do i_state=1,size_KS
       do j_state=1,size_KS
         do k_state=1,size_KS
           do l_state=1,size_KS
             num=num+1
             dset1_data(1,num,1,1,1)=Real(coul_val(i_state,j_state,k_state,&
                                          l_state))
             dset1_data(2,num,1,1,1)=Aimag(coul_val(i_state,j_state,k_state,&
                                           l_state))
           enddo
         enddo
       enddo
     enddo
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 2
     data_dims(2) = size_KS*size_KS*size_KS*size_KS
     data_dims(3) = 1
     data_dims(4) = 1
     data_dims(5) = 1
     offset(1) = 0
     offset(2) = 0
     offset(3) = q_point-1
     offset(4) = k_point-1
     offset(5) = kstrich_point-1

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_coulelement_lvl_v0_1

subroutine out_KS_vec_scalapack( eigenvec_complex_value, file_id, plist_id, &
                                 i_row, i_col, i_spin, my_k_point, dsetname)

  !  PURPOSE
  !   Output eigenvalues distributed through scalapck arrays
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  complex*16, INTENT(IN)::  eigenvec_complex_value
  integer, INTENT(IN) :: i_row
  integer, INTENT(IN) :: i_col
  integer, INTENT(IN) :: i_spin
  integer, INTENT(IN) :: my_k_point
  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  ! Write dipelement to file
  !  INPUTS
  !    o eigenvec_complex_value -- KS_eigenvalues
  !    o i_row -- i_state
  !    o i_col -- i_basis
  !    o i_spin - i_spin
  !    o my_k_point -- k-point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Data in dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(5) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     dims1 = (/1,1,1,1,2/)
     allocate(dset1_data(1,1,1,1,2))
     dset1_data(1,1,1,1,1)=REAL(eigenvec_complex_value)
     dset1_data(1,1,1,1,2)=AIMAG(eigenvec_complex_value)
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 1
     data_dims(2) = 1
     data_dims(3) = 1
     data_dims(4) = 1
     data_dims(5) = 2
     offset(1) = i_row-1
     offset(2) = i_col-1
     offset(3) = i_spin-1
     offset(4) = my_k_point-1
     offset(5) = 0

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_KS_vec_scalapack

subroutine out_ovlp_scalapack( overlap_complex_value, file_id, plist_id, &
                                 i_row, i_col, my_k_point, dsetname)

  !  PURPOSE
  !   Output overlap matrix stored in scalapack arrays, 
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  complex*16, INTENT(IN)::  overlap_complex_value !overlap_matrix(i_basis,
                                                  ! j_basis, k_point)
  integer, INTENT(IN) :: i_row  ! i_basis
  integer, INTENT(IN) :: i_col  ! j_basis
  integer, INTENT(IN) :: my_k_point ! k_point
  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  ! Write dipelement to file
  !  INPUTS
  !    o overlap_complex_value -- overlap_matrix(i_basis, j_basis, k_point)
  !    o i_row  -- i_basis
  !    o i_col  -- j_basis
  !    o my_k_point -- k_point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Data in dataset in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(4) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 4 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(4) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(4) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:)  ! Arrays to hold data

     dims1 = (/1,1,1,2/)
     allocate(dset1_data(1,1,1,2))
     dset1_data(1,1,1,1)=REAL(overlap_complex_value)
     dset1_data(1,1,1,2)=AIMAG(overlap_complex_value)
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 1
     data_dims(2) = 1
     data_dims(3) = 1
     data_dims(4) = 2
     offset(1) = i_row-1
     offset(2) = i_col-1
     offset(3) = my_k_point-1
     offset(4) = 0

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_ovlp_scalapack

subroutine out_hmlt_scalapack( hamiltonian_complex_value, file_id, plist_id, &
                                 i_row, i_col, i_spin, my_k_point, dsetname)

  !  PURPOSE
  !   Output hamiltonian matrix stored in scalapack arrays, 
  !
  ! USES
  use HDF5
  implicit none

  !  ARGUMENTS
  complex*16, INTENT(IN)::  hamiltonian_complex_value !hamiltonian_matrix(i_basis,
                                                  ! j_basis, k_point, n_spin)
  integer, INTENT(IN) :: i_row  ! i_basis
  integer, INTENT(IN) :: i_col  ! j_basis
  integer, INTENT(IN) :: i_spin  ! i_spin 
  integer, INTENT(IN) :: my_k_point ! k_point
  CHARACTER(LEN=14) :: dsetname      ! Dataset name
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  !  INPUTS
  !    o hamiltonian_complex_value -- hamiltonian_matrix(i_basis, j_basis, k_point, n_spin)
  !    o i_row  -- i_basis
  !    o i_col  -- j_basis
  !    o i_spin  -- i_spin
  !    o my_k_point -- k_point
  !    o dsetname -- Dataset name
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Data in dataset in file with file_id
  !  AUTHOR
  !      RJ Maurer, Yale University 2015
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
  !    Release version, FHI-aims (2008).
  ! SOURCE


     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier
     INTEGER(HID_T) :: memspace


     INTEGER     ::   error ! Error flag
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank


     INTEGER(HSIZE_T), DIMENSION(5) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 5 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(5) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(5) :: offset

     REAL*8, allocatable :: dset1_data(:,:,:,:,:)  ! Arrays to hold data

     dims1 = (/1,1,1,1,2/)
     allocate(dset1_data(1,1,1,1,2))
     dset1_data(1,1,1,1,1)=REAL(hamiltonian_complex_value)
     dset1_data(1,1,1,1,2)=AIMAG(hamiltonian_complex_value)
     CALL h5dopen_f(file_id, dsetname, dataset_id, error)
     data_dims(1) = 1
     data_dims(2) = 1
     data_dims(3) = 1
     data_dims(4) = 1
     data_dims(5) = 2
     offset(1) = i_row-1
     offset(2) = i_col-1
     offset(3) = my_k_point-1
     offset(4) = i_spin-1
     offset(5) = 0

     CALL h5screate_simple_f(rank,data_dims , memspace, error)
     CALL h5dget_space_f(dataset_id, dataspace_id, error)
     CALL h5sselect_hyperslab_f (dataspace_id, H5S_SELECT_SET_F, offset, &
                                 data_dims, error)
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, file_space_id = dataspace_id, &
                     mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_hmlt_scalapack

subroutine out_dimensions(file_id,plist_id)

  !  PURPOSE
  !   Write dimensions of calculation to file
  !
  ! USES
  use dimensions, only: n_basis, n_states, n_spin, n_k_points, n_atoms
  use HDF5
  implicit none

  !  ARGUMENTS
  ! 
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o dataset "dimensions" in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     INTEGER(HID_T) :: plist_id   ! Property list identifier
     CHARACTER(LEN=10) :: dsetname = "dimensions"      ! Dataset name
     INTEGER(HID_T) :: file_id        ! File identifier
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k
     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(1) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 1 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     CHARACTER(LEN=200), DIMENSION(5) :: text
     allocate(dset1_data(5))
     dims1 = (/5/)
     dset1_data(1)= n_basis
     dset1_data(2)= n_states
     dset1_data(3)= n_spin
     dset1_data(4)= n_k_points
     dset1_data(5)= n_atoms
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = 5
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
     text(1) = 'n_basis'
     text(2) = 'n_states'
     text(3) = 'n_spin'
     text(4) = 'n_k_points'
     text(5) = 'n_atoms'
     call hdf5_write_attribute(file_id, dsetname, dataset_id, 'content',&
                                data_dims(1), text)
end subroutine out_dimensions
subroutine out_basis_atom(file_id,plist_id)

  !  PURPOSE
  !   Write basis_atom to dataset
  !
  ! USES

  use basis, only: basis_atom
  use dimensions, only: n_basis
  use HDF5
  implicit none

  !  ARGUMENTS
  ! 
  !  INPUTS
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o dataset "basis_atom" in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     INTEGER(HID_T) :: plist_id   ! Property list identifier
     CHARACTER(LEN=10) :: dsetname = "basis_atom"      ! Dataset name
     INTEGER(HID_T) :: file_id        ! File identifier
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k
     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(1) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 1 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     allocate(dset1_data(n_basis))
     dims1 = (/n_basis/)
     dset1_data(1:n_basis)= basis_atom(1:n_basis)
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = n_basis
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error, xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_basis_atom

subroutine out_KS_Eigenvalue(KS_eigen, file_id,plist_id)

  !  PURPOSE
  !   Write band structure to file for all k-points 
  !
  ! USES
  use dimensions, only: n_states, n_spin, n_k_points
  use HDF5
  implicit none

  !  ARGUMENTS
  INTEGER(HID_T) :: file_id        ! File identifier
  INTEGER(HID_T) :: plist_id   ! Property list identifier
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  !  INPUTS
  !    o KS_Eigen -- KS-Eigenvalues
  !    o file_id -- File identifier
  !    o plist_id -- Property list identifier
  !  OUTPUT
  !    o Dataset "KS_Eigenvalue" in file with file_id
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
  !    Release version, FHI-aims (2008).
  ! SOURCE

     CHARACTER(LEN=13) :: dsetname = "KS_Eigenvalue"      ! Dataset name
     INTEGER(HID_T) :: dataset_id    ! Dataset1 identifier
     INTEGER(HID_T) :: dataspace_id  ! Data space identifier


     INTEGER     ::  i_k, i_spin

     INTEGER     ::   error ! Error flag
     REAL*8, allocatable :: dset1_data(:,:,:)  ! Arrays to hold data

     INTEGER(HSIZE_T), DIMENSION(3) :: dims1 ! Dataset dimensions
     INTEGER     ::   rank = 3 ! Datasets rank
     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

     allocate(dset1_data(n_states,n_spin,n_k_points))
     dims1 = (/n_states, n_spin, n_k_points/)
     do i_spin = 1, n_spin
       do i_k = 1, n_k_points
        dset1_data(1:n_states,i_spin,i_k)=KS_eigen(1:n_states, i_spin, i_k)
       enddo
     enddo
     CALL h5screate_simple_f(rank, dims1, dataspace_id, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
			      dataset_id, error)
     data_dims(1) = n_states
     data_dims(2) = n_spin     
     data_dims(3) = n_k_points
     CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset1_data, data_dims, &
                     error,xfer_prp = plist_id)
     CALL h5sclose_f(dataspace_id, error)
     CALL h5dclose_f(dataset_id, error)
     deallocate(dset1_data)
end subroutine out_KS_Eigenvalue

end module hdf5_tools
