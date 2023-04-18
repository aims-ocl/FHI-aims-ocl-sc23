#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>

#include <mpi.h>

#include <sys/types.h>
#include <sys/mman.h>

static unsigned char *ptr = 0;
static int64_t segment_size = 0;
static int my_shm_id;
static int num_shm_procs;

void aims_shm_implemented_(int *flag)
{
    *flag = 1;
}

void aims_shm_init64_(int64_t *shm_size, int *info)
{
   int shm_fd, res, i, cnt, myid, n_procs, proc_info[2], *parray;
   char shm_file_name[256];

   *info = 0;

   MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   /* Construct a unique name for the shared memory segment.
      We use actual time and PID of process 0 */

   if(myid == 0)
   {
      proc_info[0] = getpid();
      proc_info[1] = time(0);
   }
   MPI_Bcast(proc_info,2,MPI_INT,0,MPI_COMM_WORLD);

   snprintf(shm_file_name,sizeof(shm_file_name),"/fhiaims-%d-%d",proc_info[0],proc_info[1]);

   /* Open shared memory segment */

   shm_fd = shm_open(shm_file_name,O_RDWR|O_CREAT,0666);

   if(shm_fd<0) {
      fprintf(stderr,"Error opening shared memory segment %s: %s\n",shm_file_name,strerror(errno));
      *info = 1;
      return;
   }

   MPI_Barrier(MPI_COMM_WORLD);

   /* Set the requested size */

   segment_size = *shm_size;
   if(segment_size<n_procs*sizeof(int)) segment_size = n_procs*sizeof(int);

   res = ftruncate(shm_fd, segment_size);
   if(res<0) {
      fprintf(stderr,"Error setting size of shared memory segment to %.12g bytes: %s\n",
              (double)segment_size, strerror(errno));
      shm_unlink(shm_file_name);
      *info = 1;
      return;
   }

   /* Map the shared memory segment */

   ptr = (unsigned char *) mmap(0, segment_size, PROT_READ|PROT_WRITE, MAP_SHARED, shm_fd, 0);
   if(ptr == MAP_FAILED) {
      fprintf(stderr,"Error mapping shared memory segment (%.12g bytes): %s\n",
              (double)segment_size, strerror(errno));
      shm_unlink(shm_file_name);
      *info = 1;
      return;
   }

   /* Check which processes share the shared memory segment.
      We use an integer array for doing this since byte access might be dangerous */

   parray = (int*) ptr; 
   memset(parray,0,n_procs*sizeof(int));

   MPI_Barrier(MPI_COMM_WORLD);

   parray[myid] = 1;

   MPI_Barrier(MPI_COMM_WORLD);

   cnt = 0;
   for(i=0;i<n_procs;i++) {
      if(i==myid) my_shm_id = cnt;
      if(parray[i]) cnt++;
   }
   num_shm_procs = cnt;

   MPI_Barrier(MPI_COMM_WORLD);

   /* Now close the file descriptor and unlink the shared memory segment.
      It will be removed when the last process exits which is using it. */

   close(shm_fd);
   shm_unlink(shm_file_name);
}

void aims_shm_release_(void)
{
   /* Since the shared memory file is already closed and unlinked,
      the only thing to do is to unmap the memory pointer */

   if(munmap(ptr, segment_size) != 0) {
      perror("Unmap of shared memory failed");
      abort();
   }

   ptr = 0;
   segment_size = 0;

   /* Make sure all processes are ready with unmap */

   MPI_Barrier(MPI_COMM_WORLD);
}

void aims_shm_n_procs_(int *n_procs)
{
   *n_procs = num_shm_procs;
}

void aims_shm_myid_(int *myid)
{
   *myid = my_shm_id;
}

void aims_shm_get64_(void *data, int64_t *offset, int64_t *size)
{
   memcpy(data, ptr+(*offset), (size_t) *size);
}

void aims_shm_put64_(void *data, int64_t *offset, int64_t *size)
{
   memcpy(ptr+(*offset), data, (size_t) *size);
}

void aims_shm_get_ptr64_(void **c_ptr, int64_t *offset)
{
   *c_ptr = (void *) (ptr+(*offset));
}

/* 32 bit versions using plain int */

void aims_shm_init_(int *shm_size, int *info)
{
   int64_t size;

   /* If shm_size is negative, it is most likely due to an integer overflow.
      Thus we use unsigned int for converting it to 64 bits.
      This could help in situations where 2 .. 4 GB are requested */

   size = *((unsigned int*)shm_size);
   aims_shm_init64_(&size, info);
}

void aims_shm_get_(void *data, int *offset, int *size)
{
   memcpy(data, ptr+(*offset), (size_t) *size);
}

void aims_shm_put_(void *data, int *offset, int *size)
{
   memcpy(ptr+(*offset), data, (size_t) *size);
}

/* ----------------------------------------------------------------------------
**
** Below the entry points for IBM Fortran (without "_")
**
** ----------------------------------------------------------------------------
*/

void aims_shm_implemented(int *flag) {
  aims_shm_implemented_(flag);
}

void aims_shm_init64(int64_t *shm_size, int *info) {
   aims_shm_init64_(shm_size, info);
}

void aims_shm_release(void)
{
   aims_shm_release_();
}

void aims_shm_init(int *shm_size, int *info) {
   aims_shm_init_(shm_size, info);
}

void aims_shm_n_procs(int *n_procs) {
   aims_shm_n_procs_(n_procs);
}

void aims_shm_myid(int *myid) {
   aims_shm_myid_(myid);
}

void aims_shm_get64(void *data, int64_t *offset, int64_t *size) {
   aims_shm_get64_(data, offset, size);
}

void aims_shm_put64(void *data, int64_t *offset, int64_t *size) {
   aims_shm_put64_(data, offset, size);
}

void aims_shm_get_ptr64(void **c_ptr, int64_t *offset) {
   aims_shm_get_ptr64_(c_ptr, offset);
}

void aims_shm_get(void *data, int *offset, int *size) {
   aims_shm_get_(data, offset, size);
}

void aims_shm_put(void *data, int *offset, int *size) {
   aims_shm_put_(data, offset, size);
}

