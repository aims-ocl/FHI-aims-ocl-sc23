!****p* FHI-aims/aims_real
!  NAME
!    aims_real
!  SYNOPSIS
program aims_real
!  PURPOSE
!    Program to call the top-level FHI-aims subroutine.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement.
!  USE
   use mpi_utilities, only: n_tasks,myid,use_mpi,mpi_comm_global,MPI_SUCCESS,&
       initialize_mpi,finalize_mpi,aims_stop
   use localorb_io, only: localorb_allinfo
!  VARIABLES
   implicit none
   ! new MPI
   integer :: sub_tasks
   integer :: sub_class
   integer :: sub_id
   integer :: mpi_comm_subtask
   ! store current MPI
   integer :: save_myid
   integer :: save_n_tasks
   integer :: save_mpi_comm_global
   ! auxiliary
   integer :: start_job_id
   integer :: ierr
   integer :: myunit
   character*200 :: msg
   character*100 :: subfolder
   character*100, parameter :: filename = "aims.out"
!  SOURCE
   call initialize_mpi()

   ! MPI is mandatory to run multiple instances of FHI-aims
   if (.not. use_mpi) then
      call aims_stop("Multi-aims needs MPI")
   end if

   ! Read in the neccessary control information
   call read_multiaims(sub_tasks,start_job_id)

   ! Evaluate new task distribution
   sub_class = myid/sub_tasks
   sub_id = mod(myid,sub_tasks)
   write(subfolder,"(I10)") start_job_id+sub_class

   ! Give an overview of the job redistribution
   write(msg,"(3(A,I8),2A)") "Task ",myid," mapped to group ",sub_class,&
      " sub_id ",sub_id," subfolder ",trim(adjustl(subfolder))
   call localorb_allinfo(msg,unit=6)

   ! Split MPI communicator
   call MPI_Comm_split(mpi_comm_global,sub_class,sub_id,mpi_comm_subtask,ierr)
   if (ierr /= MPI_SUCCESS) then
      write(msg,"(A,I5)") "MPI_Comm_split error : ",ierr
      call aims_stop(msg,"multiaims.f90")
   end if

   call change_directory(subfolder)

   ! Save MPI meta data
   save_mpi_comm_global = mpi_comm_global
   save_myid = myid
   save_n_tasks = n_tasks

   ! Reassign myid and n_tasks in sub_class
   mpi_comm_global = mpi_comm_subtask
   myid = sub_id
   n_tasks = sub_tasks

   if (myid == 0) then
      myunit = 408+save_myid

      open(myunit,file=filename)
   end if

   call MPI_Barrier(mpi_comm_global,ierr)
   if (ierr /= MPI_SUCCESS) then
      write(msg,"(A,I5)") "MPI_Barrier error : ",ierr
      call aims_stop(msg,"multiaims.f90")
   end if

   ! Start calculation
   if (myid == 0) then
      call aims(mpi_comm_global,myunit,use_mpi)
   else
      call aims(mpi_comm_global,6,use_mpi)
   end if

   call MPI_Barrier(mpi_comm_global,ierr)
   if (ierr /= MPI_SUCCESS) then
      write(msg,"(A,I5)") "MPI_Barrier error : ",ierr
      call aims_stop(msg,"multiaims.f90")
   end if

   if (myid == 0) then
      close(myunit)
   end if

   ! Prepare for finalize
   mpi_comm_global = save_mpi_comm_global
   myid = save_myid
   n_tasks = save_n_tasks

   call MPI_Barrier(mpi_comm_global,ierr)
   if (ierr /= MPI_SUCCESS) then
      write(msg,"(A,I5)") "MPI_Barrier error : ",ierr
      call aims_stop(msg,"multiaims.f90")
   end if

   call MPI_Comm_free(mpi_comm_subtask,ierr)
   if (ierr /= MPI_SUCCESS) then
      write(msg,"(A,I5)") "MPI_comm_free error : ",ierr
      call aims_stop(msg,"multiaims.f90")
   end if

   ! Finalize mpi
   call finalize_mpi()

end program aims_real
!******
