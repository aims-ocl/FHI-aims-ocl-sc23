Module CC_cl_mpi

  Use localorb_io, only: use_unit
  Use mpi_tasks
  Use runtime_choices
  Implicit None

  Integer :: CC_mpi_comm_group
  Integer :: CC_mpi_group_size
  Integer :: CC_mpi_gid

  Integer :: CC_mpi_comm_domain
  Integer :: CC_mpi_domain_size
  Integer :: CC_mpi_did

  Double precision , dimension(:) , allocatable :: CC_mpi_real_tmp
  Double complex , dimension(:) , allocatable :: CC_mpi_cplx_tmp
  Integer , dimension(:) , allocatable :: CC_mpi_int_tmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains


  Subroutine CC_mpi_init()

  Implicit None

  Integer :: mpierr

  Integer :: color

  if (CC_n_domain.ne.1) then

    color = Mod(myid,CC_n_domain)
    Call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, myid, CC_mpi_comm_domain, mpierr)
    Call MPI_COMM_SIZE(CC_mpi_comm_domain, CC_mpi_domain_size, mpierr)
    Call MPI_COMM_RANK(CC_mpi_comm_domain, CC_mpi_did, mpierr)
  
    color = Int(myid/CC_n_domain)
    Call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, myid, CC_mpi_comm_group, mpierr)
    Call MPI_COMM_SIZE(CC_mpi_comm_group, CC_mpi_group_size, mpierr)
    Call MPI_COMM_RANK(CC_mpi_comm_group, CC_mpi_gid, mpierr)

  else

    CC_mpi_comm_domain = MPI_COMM_WORLD
    CC_mpi_domain_size = n_tasks
    CC_mpi_did = myid

    CC_mpi_comm_group = MPI_COMM_SELF
    CC_mpi_group_size = 1
    CC_mpi_gid = 0

  end if

!  print*,'CC_mpi_init',myid,MPI_COMM_WORLD,CC_mpi_did,CC_mpi_comm_domain,CC_mpi_gid,CC_mpi_comm_group

  End Subroutine CC_mpi_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_real_isend(length, vector, target_id, tag, request, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , intent(in) :: target_id
  Integer , intent(in) :: tag
  Integer , optional , intent(in) :: comm_in
  Integer , intent(out) :: request
  Double precision , dimension(length) :: vector

  Integer :: use_comm
  Integer :: ncomm, errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  if (length.gt.0) then

    Call MPI_ISEND(vector,length,MPI_DOUBLE_PRECISION, &
                      target_id,tag,use_comm,request,mpierr)

    if (mpierr.ne.0) then
      if (myid.eq.0) then
        write(use_unit,*) 'MPI err in CC_mpi_real_isend'
      end if
    end if      

  end if

  End Subroutine CC_mpi_real_isend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_real_irecv(length, vector, source_id, tag, request, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , intent(in) :: source_id
  Integer , intent(in) :: tag
  Integer , optional , intent(in) :: comm_in
  Integer , intent(out) :: request
  Double precision , dimension(length) :: vector

  Integer :: use_comm
  Integer :: ncomm, errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  if (length.gt.0) then

    Call MPI_IRECV(vector,length,MPI_DOUBLE_PRECISION, &
                      source_id,tag,use_comm,request,mpierr)

    if (mpierr.ne.0) then
      if (myid.eq.0) then
        write(use_unit,*) 'MPI err in CC_mpi_real_irecv'
      end if
    end if      

  end if

  End Subroutine CC_mpi_real_irecv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_reduce(length, vector, target_id, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , intent(in) :: target_id
  Integer , optional , intent(in) :: comm_in
  Double precision , dimension(length) , intent(inout) :: vector

  Integer(kind=8) :: cstart,cend,ccount,tmp_len,ncomm
  Integer :: use_comm
  Integer :: errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    cstart = 1
    cend = 0
  
    tmp_len = Min(length,Int(CC_MPI_max_len,8))
    Allocate(CC_mpi_real_tmp(tmp_len), stat = errnum)
    Call check_allocation(errnum,'tmp in CC_mpi_reduce')
  
    do while (cend.lt.length)
  
      ccount = length - cend
      ncomm = Min(Int(CC_MPI_max_len,8),ccount)
      cend = cend + Int(ncomm,8)
  
      Call MPI_REDUCE(vector(cstart:cend),CC_mpi_real_tmp,ncomm, &
                      MPI_DOUBLE_PRECISION,MPI_SUM,target_id,use_comm,mpierr)
      vector(cstart:cend) = CC_mpi_real_tmp(1:ncomm)
      cstart = cend + 1
  
    end do
  
    Deallocate(CC_mpi_real_tmp)

  end if

  End Subroutine CC_mpi_reduce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_allreduce(length, vector, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , optional , intent(in) :: comm_in
  Double precision , dimension(length) , intent(inout) :: vector

  Integer(kind=8) :: cstart,cend,ccount,tmp_len,ncomm
  Integer :: use_comm
  Integer :: errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    cstart = 1
    cend = 0
  
    tmp_len = Min(length,Int(CC_MPI_max_len,8))
    Allocate(CC_mpi_real_tmp(tmp_len), stat = errnum)
    Call check_allocation(errnum,'tmp in CC_mpi_allreduce')
  
    do while (cend.lt.length)
  
      ccount = length - cend
      ncomm = Min(Int(CC_MPI_max_len,8),ccount)
      cend = cend + Int(ncomm,8)
  
      Call MPI_ALLREDUCE(vector(cstart:cend),CC_mpi_real_tmp,ncomm, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,use_comm,mpierr)
  
      vector(cstart:cend) = CC_mpi_real_tmp(1:ncomm)
      cstart = cend + 1
  
    end do
  
    Deallocate(CC_mpi_real_tmp)

  end if

  End Subroutine CC_mpi_allreduce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_int_allreduce(length, vector, comm_in)

  Implicit None

  Integer , intent(in) :: length
  Integer , optional , intent(in) :: comm_in
  Integer , dimension(length) , intent(inout) :: vector

  Integer :: cstart,cend,ccount
  Integer :: use_comm
  Integer :: tmp_len, ncomm, errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    cstart = 1
    cend = 0
  
    tmp_len = Min(length,CC_MPI_max_len)
    Allocate(CC_mpi_int_tmp(tmp_len), stat = errnum)
    Call check_allocation(errnum,'tmp in CC_mpi_int_allreduce')
  
    do while (cend.lt.length)
  
      ccount = length - cend
      ncomm = Min(CC_MPI_max_len,ccount)
      cend = cend + ncomm
  
      Call MPI_ALLREDUCE(vector(cstart:cend),CC_mpi_int_tmp,ncomm, &
                         MPI_INTEGER,MPI_SUM,use_comm,mpierr)
  
      vector(cstart:cend) = CC_mpi_int_tmp(1:ncomm)
      cstart = cend + 1
  
    end do
  
    Deallocate(CC_mpi_int_tmp)

  end if

  End Subroutine CC_mpi_int_allreduce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_bcast(length, vector, source_id, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , intent(in) :: source_id
  Integer , optional , intent(in) :: comm_in
  Double precision , dimension(length) , intent(inout) :: vector

  Integer :: use_comm
  Integer :: errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    Call MPI_BCAST(vector,length,MPI_DOUBLE_PRECISION, &
                   source_id,use_comm,mpierr)

  end if

  End Subroutine CC_mpi_bcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_complex_bcast(length, vector, source_id, comm_in)

  Implicit None

  Integer(kind=8) , intent(in) :: length
  Integer , intent(in) :: source_id
  Integer , optional , intent(in) :: comm_in
  Double complex , dimension(length) , intent(inout) :: vector

  Integer (kind=8) :: cstart,cend,ccount
  Integer :: use_comm
  Integer :: ncomm, errnum, mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    Call MPI_BCAST(vector,length,MPI_DOUBLE_COMPLEX, &
                   source_id,use_comm,mpierr)

    if (mpierr.ne.0) then
      if (myid.eq.0) then
        write(use_unit,*) 'MPI err in CC_mpi_complex_bcast'
      end if
    end if

    cstart = cend + 1

  end if

  End Subroutine CC_mpi_complex_bcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_real_number(rrr,comm_in)

  Implicit None
  Integer , optional , intent(in) :: comm_in
  Double precision , intent(inout) :: rrr

  Integer :: use_comm, comm_size, mpierr
  Double precision :: r_tmp

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if (comm_size.gt.1) then

    Call MPI_ALLREDUCE(rrr,r_tmp,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,use_comm,mpierr)

    rrr = r_tmp

  end if


  End Subroutine CC_mpi_real_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_mpi_int_bcast(length, vector, source_id, comm_in)

  Implicit None

  Integer , intent(in) :: length
  Integer , intent(in) :: source_id
  Integer , optional , intent(in) :: comm_in
  Integer , dimension(length) , intent(inout) :: vector

  Integer , dimension(:) , allocatable :: tmp

  Integer :: use_comm
  Integer :: mpierr, comm_size

  if (present(comm_in)) then
    use_comm = comm_in
  else
    use_comm = MPI_COMM_WORLD
  end if

  Call MPI_COMM_SIZE(use_comm, comm_size, mpierr)

  if ((comm_size.gt.1).and.(length.gt.0)) then

    Call MPI_BCAST(vector,length,MPI_INTEGER, &
                    source_id,use_comm,mpierr)

  end if

  End Subroutine CC_mpi_int_bcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End Module CC_cl_mpi

