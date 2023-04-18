!****h* FHI-aims/distributed_hessian
!  NAME
!    distributed_hessian
!  SYNOPSIS
module distributed_hessian
  !  PURPOSE
  !    Provides variables used to distribute the Hessian matrix.
  !    Currently the distribution is a 1D block distribution - each process
  !    stores a number of rows of the whole Hessian matrix.
  !  USES
  use aims_memory_tracking, only: aims_allocate,aims_deallocate
  use dimensions, only: hess_in_file
  use localorb_io, only: use_unit,localorb_info
  use mpi_tasks
  use runtime_choices, only: relax_mode,RELAX_TRM,use_distributed_hessian
  use synchronize_mpi_basic, only: sync_vector_integer

  !  AUTHOR
  !    Victor Yu, Duke University.
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !    to the terms and conditions of the respective license agreement."
  !  HISTORY
  !    July 2018.
  !  SOURCE
  implicit none

  integer :: blk_hess ! Block size
  integer :: nrow_hess ! Number of local rows (distribtued)
  integer :: ncol_hess ! Number of local columns (not distributed)
  integer :: offset ! Defined such as global index = local index + offset
  integer :: ctxt_hess ! BLACS context
  integer :: n_worker ! Number of MPI tasks that are not idling
  integer :: comm_hess ! MPI communicator for MPI tasks that are not idling
  integer :: desc_hess(9) ! BLACS descriptor
  logical :: is_worker ! .true. if not idling

  integer, allocatable :: hess_id(:)

  integer, parameter :: password = 20180501
  character(100), parameter :: f_name = "hessian.aims"

contains
  !-----------------------------------------------------------------------------
  !****s* distributed_hessian/find_hess_dimension
  !  NAME
  !    find_hess_dimension
  !  PURPOSE
  !    Check if Hessian should be distributed, and if so, do initialization.
  !  SYNOPSIS
  subroutine find_hess_dimension(global_size)

    implicit none

    integer, intent(in) :: global_size

    integer :: id_worker
    integer :: comm_size
    integer :: ierr

    integer, external :: numroc

    character(200) :: info_str

    if (relax_mode < RELAX_TRM) then
       use_distributed_hessian = .false.
    end if

    if (.not. use_distributed_hessian) then
       blk_hess = global_size
       nrow_hess = global_size
       ncol_hess = global_size
       offset = 0
       ctxt_hess = 0
       desc_hess = 0
       is_worker = .true.
    else
       ! Find block size
       blk_hess = global_size/n_tasks

       if (blk_hess*n_tasks < global_size) then
          blk_hess = blk_hess+1
       end if

       ! Local matrix size
       nrow_hess = numroc(global_size,blk_hess,myid,0,n_tasks)
       ncol_hess = global_size

       if (nrow_hess > 0) then
          is_worker = .true.
          id_worker = 1
       else
          is_worker = .false.
          id_worker = 2
          nrow_hess = 1
       end if

       call MPI_Comm_split(mpi_comm_global,id_worker,myid,comm_hess,ierr)
       call MPI_Comm_size(comm_hess,comm_size,ierr)

       if (is_worker) then
          n_worker = comm_size
          offset = myid*blk_hess
          ctxt_hess = comm_hess

          call BLACS_Gridinit(ctxt_hess,"R",comm_size,1)
          call descinit(desc_hess,global_size,global_size,blk_hess,blk_hess,0,&
             0,ctxt_hess,max(1,nrow_hess),ierr)
       else
          n_worker = n_tasks-comm_size
          offset = 0
          ctxt_hess = 0
          desc_hess = 0
       end if
    end if ! .not. use_distributed_hessian

    if (hess_in_file .and. .not. use_distributed_hessian) then
       write(info_str,"(X,A)") &
          "*** 'hessian_file' has been found in your geometry.in, which"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") &
          "*** requires 'distributed_hessian .true.' in the control.in file."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** Please check your input files."
       call localorb_info(info_str,use_unit,"(A)")
       call aims_stop()
    end if

    if (use_distributed_hessian) then
       write(info_str,"(2X,A)") &
          "Finished initialization of distributed Hessian storage."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(2X,A,I8)") "| Global dimension:  ",ncol_hess
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(2X,A,I8)") "| BLACS block size:  ",blk_hess
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(2X,A,I8)") "| Number of workers: ",n_worker
       call localorb_info(info_str,use_unit,"(A)")
       call localorb_info("",use_unit,"(A)")
    end if

  end subroutine find_hess_dimension
  !******
  !-----------------------------------------------------------------------------
  !****s* distributed_hessian/transpose_hess
  !  NAME
  !    transpose_hess
  !  PURPOSE
  !    Transpose distributed Hessian matrix.
  !  SYNOPSIS
  subroutine transpose_hess(mat,trans_mat)

    implicit none

    real(8), intent(in)  :: mat(nrow_hess,ncol_hess)
    real(8), intent(out) :: trans_mat(nrow_hess,ncol_hess)

    integer :: data_count
    integer :: i
    integer :: ierr

    real(8), allocatable :: send_buf(:,:)
    real(8), allocatable :: recv_buf(:,:)

    if (is_worker) then
       call aims_allocate(send_buf,blk_hess,n_worker*blk_hess,"send_buf")
       call aims_allocate(recv_buf,blk_hess,n_worker*blk_hess,"recv_buf")

       send_buf = 0
       send_buf(1:nrow_hess,1:ncol_hess) = mat

       data_count = blk_hess*blk_hess

       call MPI_Alltoall(send_buf,data_count,MPI_DOUBLE_PRECISION,recv_buf,&
          data_count,MPI_DOUBLE_PRECISION,comm_hess,ierr)

       if (blk_hess > 1) then
          do i = 1,n_worker
             recv_buf(:,(i-1)*blk_hess+1:i*blk_hess) = &
                transpose(recv_buf(:,(i-1)*blk_hess+1:i*blk_hess))
          end do
       end if

       trans_mat = recv_buf(1:nrow_hess,1:ncol_hess)

       call aims_deallocate(send_buf,"send_buf")
       call aims_deallocate(recv_buf,"recv_buf")
    end if

  end subroutine transpose_hess
  !******
  !-----------------------------------------------------------------------------
  !****s* distributed_hessian/write_current_hess
  !  NAME
  !    write_current_hess
  !  SYNOPSIS
  subroutine write_current_hess(mat)
    !  PURPOSE
    !    Writes current distributed Hessian matrix to file.
    implicit none

    real(8), intent(in) :: mat(nrow_hess,ncol_hess)

    integer :: i
    integer :: i0
    integer :: i1
    integer :: f_h
    integer :: f_mode
    integer :: lsize
    integer :: ierr

    integer(8) :: f_ptr

    real(8), dimension(:), allocatable :: tmp_real

    character(200) :: info_str

    write(info_str,"(2X,3A)") &
       "Writing estimated Hessian matrix to file '",trim(f_name),"'"
    call localorb_info(info_str,use_unit,"(A)")

    if (is_worker) then
       f_mode = mpi_mode_wronly+mpi_mode_create

       call MPI_File_open(comm_hess,trim(f_name),f_mode,mpi_info_null,f_h,ierr)

       if (myid == 0) then
          f_ptr = int(0,kind=8)

          call MPI_File_write_at(f_h,f_ptr,password,1,mpi_integer4,&
             mpi_status_ignore,ierr)
       end if

       lsize = nrow_hess*ncol_hess

       call aims_allocate(tmp_real,lsize,"tmp_real")

       do i = 1,nrow_hess
          i0 = 1+(i-1)*ncol_hess
          i1 = i*ncol_hess
          tmp_real(i0:i1) = mat(i,:)
       end do

       f_ptr = int(4,kind=8)+myid*blk_hess*ncol_hess*8

       call MPI_File_write_at_all(f_h,f_ptr,tmp_real,lsize,mpi_real8,&
          mpi_status_ignore,ierr)

       call aims_deallocate(tmp_real,"tmp_real")

       call MPI_File_close(f_h,ierr)
    end if

  end subroutine write_current_hess
  !******
  !-----------------------------------------------------------------------------
  !****s* distributed_hessian/read_previous_hess
  !  NAME
  !    read_previous_hess
  !  SYNOPSIS
  subroutine read_previous_hess(mat)
    !  PURPOSE
    !    Reads previous distributed Hessian matrix from file.
    implicit none

    real(8), intent(out) :: mat(nrow_hess,ncol_hess)

    integer :: i
    integer :: i0
    integer :: i1
    integer :: f_h
    integer :: f_mode
    integer :: lsize
    integer :: tmp_int
    integer :: ierr

    integer(8) :: f_ptr

    real(8), dimension(:), allocatable :: tmp_real

    character(200) :: info_str

    write(info_str,"(2X,3A)") &
       "Reading estimated Hessian matrix from file '",trim(f_name),"'"
    call localorb_info(info_str,use_unit,"(A)")

    if (is_worker) then
       f_mode = mpi_mode_rdonly

       call MPI_File_open(comm_hess,trim(f_name),f_mode,mpi_info_null,f_h,ierr)

       if (myid == 0) then
          f_ptr = int(0,kind=8)

          call MPI_File_read_at(f_h,f_ptr,tmp_int,1,mpi_integer4,&
             mpi_status_ignore,ierr)

          if (tmp_int /= password) then
             write(info_str,"(X,3A)") "'",trim(f_name),"' is not valid."
             call localorb_info(info_str,use_unit,"(A)")
             call aims_stop()
          end if
       end if

       lsize = nrow_hess*ncol_hess

       f_ptr = int(4,kind=8)+myid*blk_hess*ncol_hess*8

       call aims_allocate(tmp_real,lsize,"tmp_real")

       call MPI_File_read_at_all(f_h,f_ptr,tmp_real,lsize,mpi_real8,&
          mpi_status_ignore,ierr)

       do i = 1,nrow_hess
          i0 = 1+(i-1)*ncol_hess
          i1 = i*ncol_hess
          mat(i,:) = tmp_real(i0:i1)
       end do

       call aims_deallocate(tmp_real,"tmp_real")

       call MPI_File_close(f_h,ierr)
    end if

  end subroutine read_previous_hess
  !******
  !-----------------------------------------------------------------------------
  !****s* distributed_hessian/check_hess_file
  !  NAME
  !    check_hess_file
  !  SYNOPSIS
  subroutine check_hess_file()
    !  PURPOSE
    !    Check the existence of "hessian.aims".
    implicit none

    logical :: f_ok
    integer :: local_ok
    integer :: bad_id
    integer :: ierr

    integer, allocatable :: global_ok(:)

    character(200) :: info_str

    write(info_str,"(2X,3A)") "Verifying the existence of the file '",&
       trim(f_name),"' on all processors"
    call localorb_info(info_str,use_unit,"(A)")

    call aims_allocate(global_ok,n_tasks,"global_ok")

    inquire(file=trim(f_name),exist=f_ok)

    if (.not. f_ok) then
       local_ok = myid
    else
       local_ok = n_tasks
    end if

    call MPI_Allgather(local_ok,1,MPI_INTEGER,global_ok,1,MPI_INTEGER,&
         mpi_comm_global,ierr)

    bad_id = minval(global_ok,1)

    if (bad_id < n_tasks) then
       if (bad_id == 0) then
          ! Not found on task 0
          write(info_str,"(X,A)") "***"
          call localorb_info(info_str,use_unit,"(A)")
          write(info_str,"(X,3A)") "*** Error - the expected file '",&
             trim(f_name),"' does not exist."
          call localorb_info(info_str,use_unit,"(A)")
       else
          ! Found on task 0 but not somewhere else
          write(info_str,"(X,A)") "***"
          call localorb_info(info_str,use_unit,"(A)")
          write(info_str,"(X,3A,I8,A)") "*** Error - the expected file '",&
             trim(f_name),"' was found on task 0 but not on task ",bad_id,&
             " and possibly other tasks."
          call localorb_info(info_str,use_unit,"(A)")
       end if

       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** This file was requested in 'geometry.in'"//&
          " by the keyword 'hessian_file'."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** Usually, this is the case when a local"//&
          " structure optimization is continued from a previous run,"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** making use of the previously accumulated"//&
          " information about the Hessian matrix."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** Having the previous Hessian matrix"//&
          " available can mean a dramatic speedup of the subsequent"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** structure optimization, which is why we"//&
          " store it, and which is why we stop the code at this"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** point to inform the user."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** There are two paths to proceed:"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,3A)") "*** (1) Copy the missing '",trim(f_name),&
          "' file from the proceding relaxation run over to the present"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***     working directory, along with the"//&
          " 'geometry.in.next_step' file used to restart the run."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** (2) Alternatively, just remove the"//&
          " 'hessian_file' keyword from the 'geometry.in.next_step' file"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***     and continue from a newly"//&
          " initialized Hessian, not from a stored version."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** Note: The Hessian is stored as a"//&
          " separate file since it will be read and stored in parallel"
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "*** throughout the following run, avoiding a"//&
          " potential memory bottleneck (per processor)."
       call localorb_info(info_str,use_unit,"(A)")
       write(info_str,"(X,A)") "***"
       call localorb_info(info_str,use_unit,"(A)")

       write(info_str,"(3A)") "The expected file '",trim(f_name),&
          "' was not found."
       call aims_stop_coll(info_str)
    end if

    call aims_deallocate(global_ok,"global_ok")

  end subroutine check_hess_file
  !******

end module distributed_hessian
!******
