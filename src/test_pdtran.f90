  subroutine test_pdtran(nsplit)

    use mpi_tasks 
    implicit none

    

    ! in/out variables

    integer, intent(in) :: nsplit


    ! local variables
    integer, parameter :: na = 663, nblk = 16
    integer            :: i

    ! variables for splitting mpi_comm_global into groups
    ! of nsplit MPI-Tasks 
    integer            :: split_row, split_col

    ! varibales for testing pdtran
    integer            :: pdtran_test_passed
    integer            :: mpierr
    integer            :: pdtran_nprocs, pdtran_myid    
    integer            :: np_cols, np_rows
    integer            :: pcol, prow
    integer            :: pdtran_comm, pdtran_blacs_ctxt
    integer            :: mxld, mxcol
    integer            :: sc_desc(9), info
    integer            :: nprow, npcol, numroc
    integer            :: istat

    real*8, allocatable :: a(:,:), z(:,:)

    split_row = myid / nsplit
    split_col = mod(myid, nsplit)
    
    call mpi_comm_split(mpi_comm_global, split_row, split_col,  &
                        pdtran_comm, mpierr)
    
    call mpi_comm_rank(pdtran_comm, pdtran_myid,  mpierr)
    call mpi_comm_size(pdtran_comm, pdtran_nprocs,mpierr)
    

    ! we want only the nsplit tasks with colour (split_row) = 0
    ! and key (split_col) = 0,1...nsplit-1 to do the test
    if (split_row .ne. 0) return
    
    do np_cols = NINT(SQRT(REAL(pdtran_nprocs))),2,-1
       if(mod(pdtran_nprocs, np_cols) == 0 ) exit
    enddo

   
   ! at the end of the above loop, pdtran_nprocs is always divisible by np_cols**2

    np_rows = pdtran_nprocs/np_cols

    pdtran_blacs_ctxt = pdtran_comm

    call BLACS_Gridinit( pdtran_blacs_ctxt, 'C', np_rows, np_cols )
    call BLACS_Gridinfo( pdtran_blacs_ctxt, nprow, npcol, prow, pcol )


   ! Determine the necessary size of the distributed matrices
   ! and set up a scalapack descriptor.

   mxld  = numroc( na, nblk, prow, 0, np_rows )
   mxcol = numroc( na, nblk, pcol, 0, np_cols )

   call descinit( sc_desc, na, na, nblk, nblk, 0, 0, pdtran_blacs_ctxt, mxld, info )


   ! if(pdtran_myid==0) then
   !    print *,' '
   !    print '(3(a,i0))', "Testing pdtran with ", pdtran_nprocs, " task(s)"
   ! endif

   allocate(a (mxld,mxcol), z (mxld,mxcol+200), stat=istat)

   if (istat .ne. 0) then
      call aims_stop("Error when allocating arrays in test_pdtran")
   endif

   a = 0
   z = 1.d99

   call pdtran(na,na,1.d0,a,1,1,sc_desc,0.d0,z,1,1,sc_desc)

   pdtran_test_passed = 0

   call BLACS_Gridexit(pdtran_blacs_ctxt)
   
   ! Check if guard value has been overwritten

   do i = mxcol+1, ubound(z,2)
      if(any(z(:,i) /= 1.d99)) pdtran_test_passed = 1
   enddo

   call MPI_ALLREDUCE(mpi_in_place, pdtran_test_passed, 1,     &
                      MPI_INTEGER, MPI_SUM, pdtran_comm, mpierr)

   if (pdtran_myid==0) then
      if (pdtran_test_passed .eq. 0) then
         ! print *, "Test passed "
         continue
      else
         write(STDERR,'(1X,A)') "* Attention. There is a known bug in some versions of Intel's math kernel library (MKL)."
         write(STDERR,'(1X,A)') "* The bug in question is in the routine 'pdtran' (parallel matrix transpose). As we have"
         write(STDERR,'(1X,A)') "* seen this bug many times, and are unable to fix it (it's Intel's, not our), we now"
         write(STDERR,'(1X,A)') "* test for buggy MKL versions and stop the code. What you should do is:"
         write(STDERR,'(1X,A)') "* EITHER recompile/run with an older MKL version (that shipped with ifort 11.0)"
         write(STDERR,'(1X,A)') "* OR compile only blacs and scalapack from scratch, on your own. "
         write(STDERR,'(1X,A)') "* The second option is harder but doable. Some instructions can be found on the aimsclub wiki."
         call aims_stop("pdtran test failed. Aborting.")
      endif
   endif


   deallocate(a, z, stat= istat)
   if (istat .ne. 0) then
 
      call aims_stop("Error when deallocating arrays in test_pdtran")
   endif
  end subroutine test_pdtran
