!****h* FHI-aims/hdf5_stubs
!  NAME
!    hdf5_stubs
!  SYNOPSIS

module hdf5 

!  PURPOSE
!  
!  Placeholders for the HDF5-parallel subroutines see http://www.hdfgroup.org
!  There are many more routines, only the ones currently used in FHI-aims have 
!  placeholders in this module.
!
!  AUTHOR
!    AIMS-TEAM
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE

  implicit none

integer, parameter :: HID_T = selected_int_kind (8)
integer, parameter :: HSIZE_T = selected_int_kind (8)
integer, parameter :: SIZE_T = selected_int_kind (8)

integer(kind=HID_T) :: H5S_SELECT_SET_F
integer(kind=HID_T) :: H5T_NATIVE_DOUBLE
integer(kind=HID_T) :: H5T_NATIVE_INTEGER
integer(kind=HID_T) :: H5T_NATIVE_CHARACTER
integer(kind=HID_T) :: H5P_FILE_ACCESS_F
integer(kind=HID_T) :: H5F_ACC_TRUNC_F
integer(kind=HID_T) :: H5F_ACC_RDWR_F
integer(kind=HID_T) :: H5P_DATASET_XFER_F
integer(kind=HID_T) :: H5FD_MPIO_INDEPENDENT_F
integer(kind=HID_T) :: H5P_DATASET_CREATE_F

interface h5dwrite_f
   module procedure h5dwrite_f_real_1d
   module procedure h5dwrite_f_real_2d
   module procedure h5dwrite_f_real_3d
   module procedure h5dwrite_f_real_4d
   module procedure h5dwrite_f_real_5d
   module procedure h5dwrite_f_real_6d
   module procedure h5dwrite_f_real_7d
   module procedure h5dwrite_f_int_1d
   module procedure h5dwrite_f_int_2d
   module procedure h5dwrite_f_int_3d
   module procedure h5dwrite_f_int_4d
   module procedure h5dwrite_f_int_5d
   module procedure h5dwrite_f_int_6d
   module procedure h5dwrite_f_int_7d
end interface h5dwrite_f

 contains
 
    subroutine no_hdf5_link_warning()
    
      use localorb_io, only: localorb_info
      use mpi_tasks, only: aims_stop
      
      character*128 :: info_str

      write(info_str,'(a)') '***** This is a call to hdf5_stubs. This means '//&
                                   'you have selected'
      call localorb_info(info_str)

      write(info_str,'(a)') '***** paralell output via hdf5 but DID NOT '//&
                                  'link against it.'
      call localorb_info(info_str)

      call aims_stop(' General STOP. You have to compile with hdf5!')
    end subroutine no_hdf5_link_warning

    SUBROUTINE h5dopen_f(loc_id, name, dset_id, hdferr) 

      IMPLICIT NONE 
      INTEGER(HID_T), INTENT(IN) :: loc_id 
      CHARACTER(LEN=*), INTENT(IN) :: name 
      INTEGER(HID_T) :: dset_id
      INTEGER :: hdferr
      
     call no_hdf5_link_warning()
      
    END SUBROUTINE h5dopen_f

    SUBROUTINE h5fopen_f (name, classtype, file_id, hdferr) 

      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: classtype
      CHARACTER(LEN=*), INTENT(IN) :: name 
      INTEGER(HID_T) :: file_id
      INTEGER :: hdferr
      
     call no_hdf5_link_warning()
      
    END SUBROUTINE h5fopen_f

    SUBROUTINE h5screate_simple_f(rank, dims, space_id, hdferr, maxdims) 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: rank 
      INTEGER(HSIZE_T), INTENT(IN) :: dims(*)
      INTEGER(HID_T) :: space_id
      INTEGER :: hdferr 
      INTEGER(HSIZE_T), OPTIONAL, INTENT(IN) :: maxdims(*)

     call no_hdf5_link_warning()
      
    END SUBROUTINE h5screate_simple_f   

    SUBROUTINE h5dget_space_f(dataset_id, dataspace_id, hdferr) 

      IMPLICIT NONE 
      INTEGER(HID_T), INTENT(IN) :: dataset_id
      INTEGER(HID_T) :: dataspace_id
      INTEGER :: hdferr
      
     call no_hdf5_link_warning()
      
    END SUBROUTINE h5dget_space_f

    SUBROUTINE h5sselect_hyperslab_f(space_id, op, start, count,&
                                     hdferr, stride, block) 
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: space_id
      INTEGER, INTENT(IN) :: op
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start 
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
      INTEGER :: hdferr
      INTEGER(HSIZE_T), DIMENSION(*), OPTIONAL, INTENT(IN) :: stride
      INTEGER(HSIZE_T), DIMENSION(*), OPTIONAL, INTENT(IN) :: block
      
     call no_hdf5_link_warning()
      
    END SUBROUTINE h5sselect_hyperslab_f

    SUBROUTINE h5dwrite_f_real_1d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_1d

    SUBROUTINE h5dwrite_f_real_2d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_2d

       SUBROUTINE h5dwrite_f_real_3d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_3d
    
       SUBROUTINE h5dwrite_f_real_4d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_4d    
    
    SUBROUTINE h5dwrite_f_real_5d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_5d
    
    SUBROUTINE h5dwrite_f_real_6d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_6d
    
    SUBROUTINE h5dwrite_f_real_7d(dset_id, mem_type_id, buf, dims, &
                               hdferr, mem_space_id, file_space_id, xfer_prp)
      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      REAL*8, INTENT(IN) :: buf(:,:,:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_real_7d
    
    SUBROUTINE h5dwrite_f_int_1d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_1d
    
    SUBROUTINE h5dwrite_f_int_2d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_2d
    
    SUBROUTINE h5dwrite_f_int_3d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_3d
    
    SUBROUTINE h5dwrite_f_int_4d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_4d

    SUBROUTINE h5dwrite_f_int_5d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_5d
    
    SUBROUTINE h5dwrite_f_int_6d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_6d
    
    SUBROUTINE h5dwrite_f_int_7d(dset_id, mem_type_id, buf, dims, & 
                              hdferr, mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE

      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER(HID_T), INTENT(IN) :: mem_type_id
      INTEGER, INTENT(IN) :: buf(:,:,:,:,:,:,:)
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims
      INTEGER(HID_T), optional, INTENT(IN) :: mem_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: file_space_id
      INTEGER(HID_T), optional, INTENT(IN) :: xfer_prp
      INTEGER, optional :: hdferr

      call no_hdf5_link_warning()
      
    END SUBROUTINE h5dwrite_f_int_7d
    
    SUBROUTINE h5sclose_f(space_id, hdferr)   
        use localorb_io
      use mpi_tasks
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: space_id
      INTEGER :: hdferr
      character*128 :: info_str
      write(info_str,'(a)') '***** This is a call to hdf5_stubs. This means '//&
                                   'you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** paralell output via hdf5 but DID NOT '//&
                                  'link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with hdf5!')
    END SUBROUTINE h5sclose_f

    SUBROUTINE h5dclose_f(dset_id, hdferr)

      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: dset_id
      INTEGER :: hdferr
      
     call no_hdf5_link_warning()
      
    END SUBROUTINE h5dclose_f

    SUBROUTINE h5dcreate_f(loc_id, name, type_id, space_id, dset_id, & 
			  hdferr, creation_prp) 
      use localorb_io
      use mpi_tasks
      IMPLICIT NONE 
      INTEGER(HID_T), INTENT(IN) :: loc_id
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER(HID_T), INTENT(IN) :: type_id
      INTEGER(HID_T), INTENT(IN) :: space_id
      INTEGER(HID_T) :: dset_id
      INTEGER :: hdferr
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: creation_prp
      
     call no_hdf5_link_warning()
    END SUBROUTINE h5dcreate_f 
                               
    SUBROUTINE h5pcreate_f(classtype, prp_id, hdferr) 
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: classtype
      INTEGER(HID_T) :: prp_id
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5pcreate_f

    SUBROUTINE h5pset_fapl_mpio_f(prp_id, comm, info, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: prp_id
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(IN) :: info
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5pset_fapl_mpio_f

    SUBROUTINE h5fcreate_f(name, access_flags, file_id, hdferr, &  
			  creation_prp, access_prp)
      use localorb_io
      use mpi_tasks
      IMPLICIT NONE 
      CHARACTER(LEN=*), INTENT(IN) :: name 
      INTEGER, INTENT(IN) :: access_flags 
      INTEGER(HID_T) :: file_id 
      INTEGER :: hdferr
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: creation_prp  
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: access_prp
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5fcreate_f

    SUBROUTINE h5pclose_f(prp_id, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: prp_id 
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5pclose_f

    SUBROUTINE h5pset_dxpl_mpio_f(prp_id, data_xfer_mode, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: prp_id
      INTEGER, INTENT(IN) :: data_xfer_mode
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5pset_dxpl_mpio_f

    SUBROUTINE h5open_f(hdferr) 
      
      IMPLICIT NONE
      INTEGER :: hdferr 
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5open_f

    SUBROUTINE h5tcopy_f(type_id, new_type_id, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: type_id
      INTEGER(HID_T) :: new_type_id 
      INTEGER :: hdferr 
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5tcopy_f

    SUBROUTINE h5tset_size_f(type_id, size, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: type_id 
      INTEGER(SIZE_T), INTENT(IN) :: size 
      INTEGER :: hdferr 
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5tset_size_f

    SUBROUTINE h5acreate_f(obj_id, name, type_id, space_id, attr_id, & 
			  hdferr, creation_prp) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: obj_id 
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER(HID_T), INTENT(IN) :: type_id 
      INTEGER(HID_T), INTENT(IN) :: space_id
      INTEGER(HID_T) :: attr_id  
      INTEGER :: hdferr      
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: creation_prp
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5acreate_f

    SUBROUTINE h5awrite_f(attr_id, memtype_id,  buf, dims, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: attr_id   
      INTEGER(HID_T), INTENT(IN) :: memtype_id
      CHARACTER, INTENT(IN) :: buf(*) 
      INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN)  :: dims 
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5awrite_f

    SUBROUTINE h5aclose_f(attr_id, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T) :: attr_id
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5aclose_f	

    SUBROUTINE h5fclose_f(file_id, hdferr) 
      
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: file_id
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5fclose_f 

    SUBROUTINE h5close_f(hdferr) 
      
      IMPLICIT NONE
      INTEGER :: hdferr
      
      call no_hdf5_link_warning()
      
    END SUBROUTINE h5close_f

    !SUBROUTINE h5pcreate_f(class, prp_id, hdferr)
    !  use localorb_io
    !  use mpi_tasks
    !  IMPLICIT NONE
    !  INTEGER(HID_T), INTENT(IN)  :: class
    !  INTEGER(HID_T), INTENT(OUT) :: prp_id
    !  INTEGER       , INTENT(OUT) :: hdferr
    !  character*128 :: info_str
    !  write(info_str,'(a)') '***** This is a call to hdf5_stubs. This means '//&
    !                               'you have selected'
    !  call localorb_info(info_str)
    !  write(info_str,'(a)') '***** paralell output via hdf5 but DID NOT '//&
    !                              'link against it.'
    !  call localorb_info(info_str)
    !  call aims_stop(' General STOP. You have to compile with hdf5!')
    !END SUBROUTINE h5pcreate_f

    SUBROUTINE h5dread_f(dset_id, mem_type_id, buf, buf_dim, hdferr, &
         mem_space_id, file_space_id, xfer_prp)

      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN) :: dset_id     ! Dataset identifier
      INTEGER(HID_T), INTENT(IN) :: mem_type_id ! Memory datatype identifier
      Integer (kind = HSIZE_T) , dimension(:) , intent(in) :: buf_dim
      Double precision , dimension(:,:,:,:) :: buf
      INTEGER :: hdferr      ! Error code
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: mem_space_id  ! Memory dataspace identfier
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: file_space_id ! File dataspace identfier
      INTEGER(HID_T), OPTIONAL, INTENT(IN) :: xfer_prp      ! Transfer property list identifier

      call no_hdf5_link_warning()

   END SUBROUTINE h5dread_f

end module hdf5
