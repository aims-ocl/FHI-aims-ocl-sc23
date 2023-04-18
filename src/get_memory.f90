subroutine get_memory()
 use mpi_tasks
 use localorb_io, only: use_unit
 implicit none

 integer, pointer, dimension(:) :: tmp
 integer, pointer, dimension(:,:) :: total_mem
 integer :: i, err, mpi_status(mpi_status_size), cpn
 character*32 :: line, in, div

 intrinsic :: sum, mod, dble, trim, adjustl

  cpn = 20 ! cores per node
           ! should at some point be replaced by an intelligent way
           ! to figure this out. for IA this should be easy by
           ! counting the numbers of hosenames that are the same,
           ! for BG I really do not now because I don't have access

  allocate(tmp(2))

  open(100,file='/proc/self/status',action='read')
  do
   read(100,'(a)',iostat=err) line
   if(err<0) exit !EOF
   if(err>0) then
     if(myid==0) then
        write(use_unit,*) 'File /proc/self/status does not exist on this architechture.'
        write(use_unit,*) 'No memory information available, returning...'
        return
     endif
   endif
   read(line,*,iostat=err) in
   select case (in)
    case('VmPeak:')
     read(line,*) in, tmp(1), div !vmpeak
    case('VmSize:')
     read(line,*) in, tmp(2), div !vmsize
   end select
  enddo
  close(100)

  allocate(total_mem(0:n_tasks-1,2))

  if(myid==0) then
   total_mem(0,:) = tmp(:)
   do i=1, n_tasks-1
    call mpi_recv(tmp, 2, mpi_integer, i, &
                  100, mpi_comm_world, mpi_status, err)
    total_mem(i,1) = tmp(1)
    total_mem(i,2) = tmp(2)
   enddo
  else
   call mpi_send(tmp, 2, mpi_integer, 0, &
                 100, mpi_comm_world, err)
  endif

  if(myid==0) then
   tmp(:) = 0d0
   write(use_unit,'(/4x,a)') 'MemStats:'
   do i=0, n_tasks-1
    write(use_unit,'(6x,a,i5.1,a,i10.3,a,i10.3,a)') &
          'Task ', i, ':  VmSize: ', total_mem(i,2), &
          ' kB, VmPeak: ', total_mem(i,1), ' kB.'
    tmp(1) = tmp(1) + total_mem(i,1)
    tmp(2) = tmp(2) + total_mem(i,2)
    if(mod(i+1,cpn)==0) then
     write(use_unit,'(8x,a,f10.3,a,f10.3,a/)') 'Usage for this node:  VmSize:', &
       dble(tmp(2))/2**10, ' MB, VmPeak: ', &
       dble(tmp(1))/2**10, ' MB'
     tmp(:) = 0d0
    endif
   enddo
   tmp(1) = sum(total_mem(:,1))
   tmp(2) = sum(total_mem(:,2))

   div = ' kB '; in = ' kB, '
   if(tmp(1)>2**10 .and. tmp(1)<=2**20) & 
     write(div,'(f10.2,a)') dble(tmp(1))/2**10, ' MB.'
   if(tmp(1)>2**20 .and. tmp(1)<=2**30) &
     write(div,'(f10.2,a)') dble(tmp(1))/2**20, ' GB.'
   if(tmp(1)>2**30) write(div,'(f10.2,a)') dble(tmp(1))/2**30, ' TB.'
   if(tmp(2)>2**10 .and. tmp(2)<=2**20) &
     write(in,'(f10.2,a)') dble(tmp(2))/2**10, ' MB, '
   if(tmp(2)>2**20) write(in,'(f10.2,a)') dble(tmp(2))/2**20,' GB, '
   write(use_unit,'(4x,a)') 'Overall total memory usage:'
   write(use_unit,'(19x,4a/)') &
           'VmSize: ', trim(adjustl(in)), ' VmPeak: ', trim(adjustl(div))
 endif
end subroutine
