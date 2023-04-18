!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium, 
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG), 
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen , 
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie, 
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaftrn, 
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition, 
!      and  
!    - IBM Deutschland GmbH
!
!
!    More information can be found here:
!    http://elpa.rzg.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the 
!    GNU Lesser General Public License as published by the Free 
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".


module ELPA1_2013

! Version 1.1.2, 2011-02-21

  implicit none

  PRIVATE ! By default, all routines contained are private

  ! The following routines are public:

  public :: get_elpa_row_col_comms_2013     ! Sets MPI row/col communicators

  public :: solve_evp_real_2013             ! Driver routine for real eigenvalue problem
  public :: solve_evp_complex_2013          ! Driver routine for complex eigenvalue problem

  public :: tridiag_real_2013               ! Transform real symmetric matrix to tridiagonal form
  public :: trans_ev_real_2013              ! Transform eigenvectors of a tridiagonal matrix back
  public :: mult_at_b_real_2013             ! Multiply real matrices A**T * B

  public :: tridiag_complex_2013            ! Transform complex hermitian matrix to tridiagonal form
  public :: trans_ev_complex_2013           ! Transform eigenvectors of a tridiagonal matrix back
  public :: mult_ah_b_complex_2013          ! Multiply complex matrices A**H * B

  public :: solve_tridi_2013                ! Solve tridiagonal eigensystem with divide and conquer method

  public :: cholesky_real_2013              ! Cholesky factorization of a real matrix
  public :: invert_trm_real_2013            ! Invert real triangular matrix

  public :: cholesky_complex_2013           ! Cholesky factorization of a complex matrix
  public :: invert_trm_complex_2013         ! Invert complex triangular matrix

  public :: local_index                ! Get local index of a block cyclic distributed matrix
  public :: least_common_multiple      ! Get least common multiple

  public :: hh_transform_real
  public :: hh_transform_complex

!-------------------------------------------------------------------------------

  ! Timing results, set by every call to solve_evp_xxx

  real*8, public :: time_evp_fwd    ! forward transformations (to tridiagonal form)
  real*8, public :: time_evp_solve  ! time for solving the tridiagonal system
  real*8, public :: time_evp_back   ! time for back transformations of eigenvectors

  ! Set elpa_print_times to .true. for explicit timing outputs

  logical, public :: elpa_print_times = .false.

!-------------------------------------------------------------------------------

  include 'mpif.h'

contains

!-------------------------------------------------------------------------------

subroutine get_elpa_row_col_comms_2013(mpi_comm_global, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
! get_elpa_row_col_comms:
! All ELPA routines need MPI communicators for communicating within
! rows or columns of processes, these are set here.
! mpi_comm_rows/mpi_comm_cols can be free'd with MPI_Comm_free if not used any more.
!
!  Parameters
!
!  mpi_comm_global   Global communicator for the calculations (in)
!
!  my_prow           Row coordinate of the calling process in the process grid (in)
!
!  my_pcol           Column coordinate of the calling process in the process grid (in)
!
!  mpi_comm_rows     Communicator for communicating within rows of processes (out)
!
!  mpi_comm_cols     Communicator for communicating within columns of processes (out)
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: mpi_comm_global, my_prow, my_pcol
   integer, intent(out) :: mpi_comm_rows, mpi_comm_cols

   integer :: mpierr

   ! mpi_comm_rows is used for communicating WITHIN rows, i.e. all processes
   ! having the same column coordinate share one mpi_comm_rows.
   ! So the "color" for splitting is my_pcol and the "key" is my row coordinate.
   ! Analogous for mpi_comm_cols

   call mpi_comm_split(mpi_comm_global,my_pcol,my_prow,mpi_comm_rows,mpierr)
   call mpi_comm_split(mpi_comm_global,my_prow,my_pcol,mpi_comm_cols,mpierr)

end subroutine get_elpa_row_col_comms_2013

!-------------------------------------------------------------------------------

subroutine solve_evp_real_2013(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_global)

!-------------------------------------------------------------------------------
!  solve_evp_real: Solves the real eigenvalue problem
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed.
!              The smallest nev eigenvalues/eigenvectors are calculated.
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_global optional and not used here, just to be consistent with the
!                  2017 ELPA Version
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 :: a(lda,*), ev(na), q(ldq,*)

   integer my_prow, my_pcol, mpierr
   real*8, allocatable :: e(:), tau(:)
   real*8 ttt0, ttt1
   integer, optional :: mpi_comm_global

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)

   allocate(e(na), tau(na))

   ttt0 = MPI_Wtime()
   call tridiag_real_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time tridiag_real :',ttt1-ttt0
   time_evp_fwd = ttt1-ttt0

   ttt0 = MPI_Wtime()
   call solve_tridi_2013(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time solve_tridi  :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0

   ttt0 = MPI_Wtime()
   call trans_ev_real_2013(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time trans_ev_real:',ttt1-ttt0
   time_evp_back = ttt1-ttt0

   deallocate(e, tau)

end subroutine solve_evp_real_2013

!-------------------------------------------------------------------------------


subroutine solve_evp_complex_2013(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  solve_evp_complex: Solves the complex eigenvalue problem
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!              The smallest nev eigenvalues/eigenvectors are calculated.
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   complex*16 :: a(lda,*), q(ldq,*)
   real*8 :: ev(na)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_rows, l_cols, l_cols_nev
   real*8, allocatable :: q_real(:,:), e(:)
   complex*16, allocatable :: tau(:)
   real*8 ttt0, ttt1

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

   l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

   allocate(e(na), tau(na))
   allocate(q_real(l_rows,l_cols))

   ttt0 = MPI_Wtime()
   call tridiag_complex_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time tridiag_complex :',ttt1-ttt0
   time_evp_fwd = ttt1-ttt0

   ttt0 = MPI_Wtime()
   call solve_tridi_2013(na, nev, ev, e, q_real, l_rows, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time solve_tridi     :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0

   ttt0 = MPI_Wtime()
   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   call trans_ev_complex_2013(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) print *,'Time trans_ev_complex:',ttt1-ttt0
   time_evp_back = ttt1-ttt0

   deallocate(q_real)
   deallocate(e, tau)

end subroutine solve_evp_complex_2013

!-------------------------------------------------------------------------------

subroutine tridiag_real_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, d, e, tau)

!-------------------------------------------------------------------------------
!  tridiag_real: Reduces a distributed symmetric matrix to tridiagonal form
!                (like Scalapack Routine PDSYTRD)
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the Householder vectors
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  d(na)       Diagonal elements (returned), identical on all processors
!
!  e(na)       Off-Diagonal elements (returned), identical on all processors
!
!  tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*), d(na), e(na), tau(na)

   integer, parameter :: max_stored_rows = 32

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, nstor
   integer istep, i, j, lcs, lce, lrs, lre
   integer tile_size, l_rows_tile, l_cols_tile

   integer my_thread, n_threads, max_threads, n_iter
   integer omp_get_thread_num, omp_get_num_threads, omp_get_max_threads

   real*8 vav, vnorm2, x, aux(2*max_stored_rows), aux1(2), aux2(2), vrl, xf

   real*8, allocatable:: tmp(:), vr(:), vc(:), ur(:), uc(:), vur(:,:), uvc(:,:)

   real*8, allocatable:: ur_p(:,:), uc_p(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = (totalblocks-1)/np_cols + 1

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   allocate(tmp(MAX(max_local_rows,max_local_cols)))
   allocate(vr(max_local_rows+1))
   allocate(ur(max_local_rows))
   allocate(vc(max_local_cols))
   allocate(uc(max_local_cols))

   max_threads = omp_get_max_threads()

   allocate(ur_p(max_local_rows,0:max_threads-1))
   allocate(uc_p(max_local_cols,0:max_threads-1))


   tmp = 0
   vr = 0
   ur = 0
   vc = 0
   uc = 0

   allocate(vur(max_local_rows,2*max_stored_rows))
   allocate(uvc(max_local_cols,2*max_stored_rows))

   d(:) = 0
   e(:) = 0
   tau(:) = 0

   nstor = 0

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a
   if(my_prow==prow(na) .and. my_pcol==pcol(na)) d(na) = a(l_rows,l_cols)

   do istep=na,3,-1

      ! Calculate number of local rows and columns of the still remaining matrix
      ! on the local processor

      l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
      l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

      ! Calculate vector for Householder transformation on all procs
      ! owning column istep

      if(my_pcol==pcol(istep)) then

         ! Get vector to be transformed; distribute last element and norm of
         ! remaining elements to all procs in current column

         vr(1:l_rows) = a(1:l_rows,l_cols+1)
         if(nstor>0 .and. l_rows>0) then
            call DGEMV('N',l_rows,2*nstor,1.d0,vur,ubound(vur,1), &
                       uvc(l_cols+1,1),ubound(uvc,1),1.d0,vr,1)
         endif

         if(my_prow==prow(istep-1)) then
            aux1(1) = dot_product(vr(1:l_rows-1),vr(1:l_rows-1))
            aux1(2) = vr(l_rows)
         else
            aux1(1) = dot_product(vr(1:l_rows),vr(1:l_rows))
            aux1(2) = 0.
         endif

         call mpi_allreduce(aux1,aux2,2,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

         vnorm2 = aux2(1)
         vrl    = aux2(2)

         ! Householder transformation

         call hh_transform_real(vrl, vnorm2, xf, tau(istep))

         ! Scale vr and store Householder vector for back transformation

         vr(1:l_rows) = vr(1:l_rows) * xf
         if(my_prow==prow(istep-1)) then
            vr(l_rows) = 1.
            e(istep-1) = vrl
         endif
         a(1:l_rows,l_cols+1) = vr(1:l_rows) ! store Householder vector for back transformation

      endif

      ! Broadcast the Householder vector (and tau) along columns

      if(my_pcol==pcol(istep)) vr(l_rows+1) = tau(istep)
      call MPI_Bcast(vr,l_rows+1,MPI_REAL8,pcol(istep),mpi_comm_cols,mpierr)
      tau(istep) =  vr(l_rows+1)

      ! Transpose Householder vector vr -> vc

      call elpa_transpose_vectors  (vr, ubound(vr,1), mpi_comm_rows, &
                                    vc, ubound(vc,1), mpi_comm_cols, &
                                    1, istep-1, 1, nblk)


      ! Calculate u = (A + VU**T + UV**T)*v

      ! For cache efficiency, we use only the upper half of the matrix tiles for this,
      ! thus the result is partly in uc(:) and partly in ur(:)

      uc(1:l_cols) = 0
      ur(1:l_rows) = 0
      if(l_rows>0 .and. l_cols>0) then

!$OMP PARALLEL PRIVATE(my_thread,n_threads,n_iter,i,lcs,lce,j,lrs,lre)

         my_thread = omp_get_thread_num()
         n_threads = omp_get_num_threads()

         n_iter = 0

         uc_p(1:l_cols,my_thread) = 0.
         ur_p(1:l_rows,my_thread) = 0.

         do i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if(lce<lcs) cycle
            do j=0,i
               lrs = j*l_rows_tile+1
               lre = min(l_rows,(j+1)*l_rows_tile)
               if(lre<lrs) cycle

               if(mod(n_iter,n_threads) == my_thread) then
                 call DGEMV('T',lre-lrs+1,lce-lcs+1,1.d0,a(lrs,lcs),lda,vr(lrs),1,1.d0,uc_p(lcs,my_thread),1)
                 if(i/=j) call DGEMV('N',lre-lrs+1,lce-lcs+1,1.d0,a(lrs,lcs),lda,vc(lcs),1,1.d0,ur_p(lrs,my_thread),1)
               endif
               n_iter = n_iter+1

            enddo
         enddo

!$OMP END PARALLEL
      
         do i=0,max_threads-1
            uc(1:l_cols) = uc(1:l_cols) + uc_p(1:l_cols,i)
            ur(1:l_rows) = ur(1:l_rows) + ur_p(1:l_rows,i)
         enddo

         if(nstor>0) then
            call DGEMV('T',l_rows,2*nstor,1.d0,vur,ubound(vur,1),vr,1,0.d0,aux,1)
            call DGEMV('N',l_cols,2*nstor,1.d0,uvc,ubound(uvc,1),aux,1,1.d0,uc,1)
         endif

      endif

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      if(tile_size < istep-1) then
         call elpa_reduce_add_vectors  (ur, ubound(ur,1), mpi_comm_rows, &
                                        uc, ubound(uc,1), mpi_comm_cols, &
                                        istep-1, 1, nblk)
      endif

      ! Sum up all the uc(:) parts, transpose uc -> ur

      if(l_cols>0) then
         tmp(1:l_cols) = uc(1:l_cols)
         call mpi_allreduce(tmp,uc,l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
      endif

      call elpa_transpose_vectors  (uc, ubound(uc,1), mpi_comm_cols, &
                                    ur, ubound(ur,1), mpi_comm_rows, &
                                    1, istep-1, 1, nblk)

      ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

      x = 0
      if(l_cols>0) x = dot_product(vc(1:l_cols),uc(1:l_cols))
      call mpi_allreduce(x,vav,1,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)

      ! store u and v in the matrices U and V
      ! these matrices are stored combined in one here

      do j=1,l_rows
         vur(j,2*nstor+1) = tau(istep)*vr(j)
         vur(j,2*nstor+2) = 0.5*tau(istep)*vav*vr(j) - ur(j)
      enddo
      do j=1,l_cols
         uvc(j,2*nstor+1) = 0.5*tau(istep)*vav*vc(j) - uc(j)
         uvc(j,2*nstor+2) = tau(istep)*vc(j)
      enddo

      nstor = nstor+1

      ! If the limit of max_stored_rows is reached, calculate A + VU**T + UV**T

      if(nstor==max_stored_rows .or. istep==3) then

         do i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = 1
            lre = min(l_rows,(i+1)*l_rows_tile)
            if(lce<lcs .or. lre<lrs) cycle
            call dgemm('N','T',lre-lrs+1,lce-lcs+1,2*nstor,1.d0, &
                       vur(lrs,1),ubound(vur,1),uvc(lcs,1),ubound(uvc,1), &
                       1.d0,a(lrs,lcs),lda)
         enddo

         nstor = 0

      endif

      if(my_prow==prow(istep-1) .and. my_pcol==pcol(istep-1)) then
         if(nstor>0) a(l_rows,l_cols) = a(l_rows,l_cols) &
                        + dot_product(vur(l_rows,1:2*nstor),uvc(l_cols,1:2*nstor))
         d(istep-1) = a(l_rows,l_cols)
      endif

   enddo

   ! Store e(1) and d(1)

   if(my_prow==prow(1) .and. my_pcol==pcol(2)) e(1) = a(1,l_cols) ! use last l_cols value of loop above
   if(my_prow==prow(1) .and. my_pcol==pcol(1)) d(1) = a(1,1)

   deallocate(tmp, vr, ur, vc, uc, vur, uvc)

   ! distribute the arrays d and e to all processors

   allocate(tmp(na))
   tmp = d
   call mpi_allreduce(tmp,d,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmp = d
   call mpi_allreduce(tmp,d,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   tmp = e
   call mpi_allreduce(tmp,e,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmp = e
   call mpi_allreduce(tmp,e,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   deallocate(tmp)

end subroutine tridiag_real_2013

!-------------------------------------------------------------------------------

subroutine trans_ev_real_2013(na, nqc, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_real: Transforms the eigenvectors of a tridiagonal matrix back
!                 to the eigenvectors of the original matrix
!                 (like Scalapack Routine PDORMTR)
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tau(na)     Factors of the Householder vectors
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, nqc, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*), q(ldq,*), tau(na)

   integer :: max_stored_rows

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, nstor
   integer istep, i, n, nc, ic, ics, ice, nb, cur_pcol

   real*8, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)
   real*8, allocatable:: tmat(:,:), h1(:), h2(:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk


   max_stored_rows = (63/nblk+1)*nblk

   allocate(tmat(max_stored_rows,max_stored_rows))
   allocate(h1(max_stored_rows*max_stored_rows))
   allocate(h2(max_stored_rows*max_stored_rows))
   allocate(tmp1(max_local_cols*max_stored_rows))
   allocate(tmp2(max_local_cols*max_stored_rows))
   allocate(hvb(max_local_rows*nblk))
   allocate(hvm(max_local_rows,max_stored_rows))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   nstor = 0

   do istep=1,na,nblk

      ics = MAX(istep,3)
      ice = MIN(istep+nblk-1,na)
      if(ice<ics) cycle

      cur_pcol = pcol(istep)

      nb = 0
      do ic=ics,ice

         l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder vector
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector


         if(my_pcol==cur_pcol) then
            hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)
            if(my_prow==prow(ic-1)) then
               hvb(nb+l_rows) = 1.
            endif
         endif

         nb = nb+l_rows
      enddo

      if(nb>0) &
         call MPI_Bcast(hvb,nb,MPI_REAL8,cur_pcol,mpi_comm_cols,mpierr)

      nb = 0
      do ic=ics,ice
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
         hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
         nstor = nstor+1
         nb = nb+l_rows
      enddo

      ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
      if(nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32)) then

         ! Calculate scalar products of stored vectors.
         ! This can be done in different ways, we use dsyrk

         tmat = 0
         if(l_rows>0) &
            call dsyrk('U','T',nstor,l_rows,1.d0,hvm,ubound(hvm,1),0.d0,tmat,max_stored_rows)

         nc = 0
         do n=1,nstor-1
            h1(nc+1:nc+n) = tmat(1:n,n+1)
            nc = nc+n
         enddo

         if(nc>0) call mpi_allreduce(h1,h2,nc,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)

         ! Calculate triangular matrix T

         nc = 0
         tmat(1,1) = tau(ice-nstor+1)
         do n=1,nstor-1
            call dtrmv('L','T','N',n,tmat,max_stored_rows,h2(nc+1),1)
            tmat(n+1,1:n) = -h2(nc+1:nc+n)*tau(ice-nstor+n+1)
            tmat(n+1,n+1) = tau(ice-nstor+n+1)
            nc = nc+n
         enddo

         ! Q = Q - V * T * V**T * Q

         if(l_rows>0) then
            call dgemm('T','N',nstor,l_cols,l_rows,1.d0,hvm,ubound(hvm,1), &
                       q,ldq,0.d0,tmp1,nstor)
         else
            tmp1(1:l_cols*nstor) = 0
         endif
         call mpi_allreduce(tmp1,tmp2,nstor*l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
         if(l_rows>0) then
            call dtrmm('L','L','N','N',nstor,l_cols,1.0d0,tmat,max_stored_rows,tmp2,nstor)
            call dgemm('N','N',l_rows,l_cols,nstor,-1.d0,hvm,ubound(hvm,1), &
                       tmp2,nstor,1.d0,q,ldq)
         endif
         nstor = 0
      endif

   enddo

   deallocate(tmat, h1, h2, tmp1, tmp2, hvb, hvm)


end subroutine trans_ev_real_2013

!-------------------------------------------------------------------------------

subroutine mult_at_b_real_2013(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

!-------------------------------------------------------------------------------
!  mult_at_b_real:  Performs C := A**T * B
!
!      where:  A is a square matrix (na,na) which is optionally upper or lower triangular
!              B is a (na,ncb) matrix
!              C is a (na,ncb) matrix where optionally only the upper or lower
!              triangle may be computed
!
!  Parameters
!
!  uplo_a      'U' if A is upper triangular
!              'L' if A is lower triangular
!              anything else if A is a full matrix
!              Please note: This pertains to the original A (as set in the calling program)
!              whereas the transpose of A is used for calculations
!              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!              i.e. it may contain arbitrary numbers
!
!  uplo_c      'U' if only the upper diagonal part of C is needed
!              'L' if only the upper diagonal part of C is needed
!              anything else if the full matrix C is needed
!              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!              written to a certain extent, i.e. one shouldn't rely on the content there!
!
!  na          Number of rows/columns of A, number of rows of B and C
!
!  ncb         Number of columns  of B and C
!
!  a           Matrix A
!
!  lda         Leading dimension of a
!
!  b           Matrix B
!
!  ldb         Leading dimension of b
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  c           Matrix C
!
!  ldc         Leading dimension of c
!
!-------------------------------------------------------------------------------

   implicit none

   character*1 uplo_a, uplo_c

   integer na, ncb, lda, ldb, nblk, mpi_comm_rows, mpi_comm_cols, ldc
   real*8 a(lda,*), b(ldb,*), c(ldc,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_rows_np
   integer np, n, nb, nblk_mult, lrs, lre, lcs, lce
   integer gcol_min, gcol, goff
   integer nstor, nr_done, noff, np_bc, n_aux_bc, nvals
   integer, allocatable :: lrs_save(:), lre_save(:)

   logical a_lower, a_upper, c_lower, c_upper

   real*8, allocatable:: aux_mat(:,:), aux_bc(:), tmp1(:,:), tmp2(:,:)


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
   l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

   ! Block factor for matrix multiplications, must be a multiple of nblk

   if(na/np_rows<=256) then
      nblk_mult = (31/nblk+1)*nblk
   else
      nblk_mult = (63/nblk+1)*nblk
   endif

   allocate(aux_mat(l_rows,nblk_mult))
   allocate(aux_bc(l_rows*nblk))
   allocate(lrs_save(nblk))
   allocate(lre_save(nblk))

   a_lower = .false.
   a_upper = .false.
   c_lower = .false.
   c_upper = .false.

   if(uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
   if(uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
   if(uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
   if(uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

   ! Build up the result matrix by processor rows

   do np = 0, np_rows-1

      ! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0   ! Number of columns stored in aux_mat

      ! Loop over the blocks on row np

      do nb=0,(l_rows_np-1)/nblk

         goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

         ! Get the processor column which owns this block (A is transposed, so we need the column)
         ! and the offset in blocks within this column.
         ! The corresponding block column in A is then broadcast to all for multiplication with B

         np_bc = MOD(goff,np_cols)
         noff = goff/np_cols
         n_aux_bc = 0

         ! Gather up the complete block column of A on the owner

         do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

            gcol = goff*nblk + n ! global column corresponding to n
            if(nstor==0 .and. n==1) gcol_min = gcol

            lrs = 1       ! 1st local row number for broadcast
            lre = l_rows  ! last local row number for broadcast
            if(a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
            if(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            if(lrs<=lre) then
               nvals = lre-lrs+1
               if(my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
               n_aux_bc = n_aux_bc + nvals
            endif

            lrs_save(n) = lrs
            lre_save(n) = lre

         enddo

         ! Broadcast block column

         call MPI_Bcast(aux_bc,n_aux_bc,MPI_REAL8,np_bc,mpi_comm_cols,mpierr)

         ! Insert what we got in aux_mat

         n_aux_bc = 0
         do n = 1, min(l_rows_np-nb*nblk,nblk)
            nstor = nstor+1
            lrs = lrs_save(n)
            lre = lre_save(n)
            if(lrs<=lre) then
               nvals = lre-lrs+1
               aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
               n_aux_bc = n_aux_bc + nvals
            endif
         enddo

         ! If we got nblk_mult columns in aux_mat or this is the last block
         ! do the matrix multiplication

         if(nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

            lrs = 1       ! 1st local row number for multiply
            lre = l_rows  ! last local row number for multiply
            if(a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
            if(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            lcs = 1       ! 1st local col number for multiply
            lce = l_cols  ! last local col number for multiply
            if(c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
            if(c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

            if(lcs<=lce) then
               allocate(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce))
               if(lrs<=lre) then
                  call dgemm('T','N',nstor,lce-lcs+1,lre-lrs+1,1.d0,aux_mat(lrs,1),ubound(aux_mat,1), &
                             b(lrs,lcs),ldb,0.d0,tmp1,nstor)
               else
                  tmp1 = 0
               endif

               ! Sum up the results and send to processor row np
               call mpi_reduce(tmp1,tmp2,nstor*(lce-lcs+1),MPI_REAL8,MPI_SUM,np,mpi_comm_rows,mpierr)

               ! Put the result into C
               if(my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

               deallocate(tmp1,tmp2)
            endif

            nr_done = nr_done+nstor
            nstor=0
            aux_mat(:,:)=0
         endif
      enddo
   enddo

   deallocate(aux_mat, aux_bc, lrs_save, lre_save)

end subroutine mult_at_b_real_2013

!-------------------------------------------------------------------------------

subroutine tridiag_complex_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, d, e, tau)

!-------------------------------------------------------------------------------
!  tridiag_complex: Reduces a distributed hermitian matrix to tridiagonal form
!                   (like Scalapack Routine PZHETRD)
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to PZHETRD, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the Householder vectors
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  d(na)       Diagonal elements (returned), identical on all processors
!
!  e(na)       Off-Diagonal elements (returned), identical on all processors
!
!  tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), tau(na)
   real*8 d(na), e(na)

   integer, parameter :: max_stored_rows = 32

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, nstor
   integer istep, i, j, lcs, lce, lrs, lre
   integer tile_size, l_rows_tile, l_cols_tile

   integer my_thread, n_threads, max_threads, n_iter
   integer omp_get_thread_num, omp_get_num_threads, omp_get_max_threads


   real*8 vnorm2
   complex*16 vav, xc, aux(2*max_stored_rows),  aux1(2), aux2(2), vrl, xf

   complex*16, allocatable:: tmp(:), vr(:), vc(:), ur(:), uc(:), vur(:,:), uvc(:,:)

   complex*16, allocatable:: ur_p(:,:), uc_p(:,:)

   real*8, allocatable:: tmpr(:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = (totalblocks-1)/np_cols + 1

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   allocate(tmp(MAX(max_local_rows,max_local_cols)))
   allocate(vr(max_local_rows+1))
   allocate(ur(max_local_rows))
   allocate(vc(max_local_cols))
   allocate(uc(max_local_cols))

   max_threads = omp_get_max_threads()

   allocate(ur_p(max_local_rows,0:max_threads-1))
   allocate(uc_p(max_local_cols,0:max_threads-1))


   tmp = 0
   vr = 0
   ur = 0
   vc = 0
   uc = 0

   allocate(vur(max_local_rows,2*max_stored_rows))
   allocate(uvc(max_local_cols,2*max_stored_rows))

   d(:) = 0
   e(:) = 0
   tau(:) = 0

   nstor = 0

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a
   if(my_prow==prow(na) .and. my_pcol==pcol(na)) d(na) = a(l_rows,l_cols)

   do istep=na,3,-1

      ! Calculate number of local rows and columns of the still remaining matrix
      ! on the local processor

      l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
      l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

      ! Calculate vector for Householder transformation on all procs
      ! owning column istep

      if(my_pcol==pcol(istep)) then

         ! Get vector to be transformed; distribute last element and norm of
         ! remaining elements to all procs in current column

         vr(1:l_rows) = a(1:l_rows,l_cols+1)
         if(nstor>0 .and. l_rows>0) then
            aux(1:2*nstor) = conjg(uvc(l_cols+1,1:2*nstor))
            call ZGEMV('N',l_rows,2*nstor,CONE,vur,ubound(vur,1), &
                       aux,1,CONE,vr,1)
         endif

         if(my_prow==prow(istep-1)) then
            aux1(1) = dot_product(vr(1:l_rows-1),vr(1:l_rows-1))
            aux1(2) = vr(l_rows)
         else
            aux1(1) = dot_product(vr(1:l_rows),vr(1:l_rows))
            aux1(2) = 0.
         endif

         call mpi_allreduce(aux1,aux2,2,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         vnorm2 = aux2(1)
         vrl    = aux2(2)

         ! Householder transformation

         call hh_transform_complex(vrl, vnorm2, xf, tau(istep))

         ! Scale vr and store Householder vector for back transformation

         vr(1:l_rows) = vr(1:l_rows) * xf
         if(my_prow==prow(istep-1)) then
            vr(l_rows) = 1.
            e(istep-1) = vrl
         endif
         a(1:l_rows,l_cols+1) = vr(1:l_rows) ! store Householder vector for back transformation

      endif

      ! Broadcast the Householder vector (and tau) along columns

      if(my_pcol==pcol(istep)) vr(l_rows+1) = tau(istep)
      call MPI_Bcast(vr,l_rows+1,MPI_DOUBLE_COMPLEX,pcol(istep),mpi_comm_cols,mpierr)
      tau(istep) =  vr(l_rows+1)

      ! Transpose Householder vector vr -> vc

      call elpa_transpose_vectors  (vr, 2*ubound(vr,1), mpi_comm_rows, &
                                    vc, 2*ubound(vc,1), mpi_comm_cols, &
                                    1, 2*(istep-1), 1, 2*nblk)

      ! Calculate u = (A + VU**T + UV**T)*v

      ! For cache efficiency, we use only the upper half of the matrix tiles for this,
      ! thus the result is partly in uc(:) and partly in ur(:)

      uc(1:l_cols) = 0
      ur(1:l_rows) = 0
      if(l_rows>0 .and. l_cols>0) then

!$OMP PARALLEL PRIVATE(my_thread,n_threads,n_iter,i,lcs,lce,j,lrs,lre)

         my_thread = omp_get_thread_num()
         n_threads = omp_get_num_threads()

         n_iter = 0

         uc_p(1:l_cols,my_thread) = 0.
         ur_p(1:l_rows,my_thread) = 0.

         do i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if(lce<lcs) cycle
            do j=0,i
               lrs = j*l_rows_tile+1
               lre = min(l_rows,(j+1)*l_rows_tile)
               if(lre<lrs) cycle

               if(mod(n_iter,n_threads) == my_thread) then
                  call ZGEMV('C',lre-lrs+1,lce-lcs+1,CONE,a(lrs,lcs),lda,vr(lrs),1,CONE,uc_p(lcs,my_thread),1)
                  if(i/=j) call ZGEMV('N',lre-lrs+1,lce-lcs+1,CONE,a(lrs,lcs),lda,vc(lcs),1,CONE,ur_p(lrs,my_thread),1)
               endif
               n_iter = n_iter+1

            enddo
         enddo

!$OMP END PARALLEL

         do i=0,max_threads-1
            uc(1:l_cols) = uc(1:l_cols) + uc_p(1:l_cols,i)
            ur(1:l_rows) = ur(1:l_rows) + ur_p(1:l_rows,i)
         enddo


         if(nstor>0) then
            call ZGEMV('C',l_rows,2*nstor,CONE,vur,ubound(vur,1),vr,1,CZERO,aux,1)
            call ZGEMV('N',l_cols,2*nstor,CONE,uvc,ubound(uvc,1),aux,1,CONE,uc,1)
         endif

      endif

      ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      if(tile_size < istep-1) then
         call elpa_reduce_add_vectors  (ur, 2*ubound(ur,1), mpi_comm_rows, &
                                        uc, 2*ubound(uc,1), mpi_comm_cols, &
                                        2*(istep-1), 1, 2*nblk)
      endif

      ! Sum up all the uc(:) parts, transpose uc -> ur

      if(l_cols>0) then
         tmp(1:l_cols) = uc(1:l_cols)
         call mpi_allreduce(tmp,uc,l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
      endif

      call elpa_transpose_vectors  (uc, 2*ubound(uc,1), mpi_comm_cols, &
                                    ur, 2*ubound(ur,1), mpi_comm_rows, &
                                    1, 2*(istep-1), 1, 2*nblk)

      ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

      xc = 0
      if(l_cols>0) xc = dot_product(vc(1:l_cols),uc(1:l_cols))
      call mpi_allreduce(xc,vav,1,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_cols,mpierr)

      ! store u and v in the matrices U and V
      ! these matrices are stored combined in one here

      do j=1,l_rows
         vur(j,2*nstor+1) = conjg(tau(istep))*vr(j)
         vur(j,2*nstor+2) = 0.5*conjg(tau(istep))*vav*vr(j) - ur(j)
      enddo
      do j=1,l_cols
         uvc(j,2*nstor+1) = 0.5*conjg(tau(istep))*vav*vc(j) - uc(j)
         uvc(j,2*nstor+2) = conjg(tau(istep))*vc(j)
      enddo

      nstor = nstor+1

      ! If the limit of max_stored_rows is reached, calculate A + VU**T + UV**T

      if(nstor==max_stored_rows .or. istep==3) then

         do i=0,(istep-2)/tile_size
            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = 1
            lre = min(l_rows,(i+1)*l_rows_tile)
            if(lce<lcs .or. lre<lrs) cycle
            call ZGEMM('N','C',lre-lrs+1,lce-lcs+1,2*nstor,CONE, &
                       vur(lrs,1),ubound(vur,1),uvc(lcs,1),ubound(uvc,1), &
                       CONE,a(lrs,lcs),lda)
         enddo

         nstor = 0

      endif

      if(my_prow==prow(istep-1) .and. my_pcol==pcol(istep-1)) then
         if(nstor>0) a(l_rows,l_cols) = a(l_rows,l_cols) &
                        + dot_product(vur(l_rows,1:2*nstor),uvc(l_cols,1:2*nstor))
         d(istep-1) = a(l_rows,l_cols)
      endif

   enddo

   ! Store e(1) and d(1)

   if(my_pcol==pcol(2)) then
      if(my_prow==prow(1)) then
         ! We use last l_cols value of loop above
         vrl = a(1,l_cols)
         call hh_transform_complex(vrl, 0.d0, xf, tau(2))
         e(1) = vrl
         a(1,l_cols) = 1. ! for consistency only
      endif
      call mpi_bcast(tau(2),1,MPI_DOUBLE_COMPLEX,prow(1),mpi_comm_rows,mpierr)
   endif
   call mpi_bcast(tau(2),1,MPI_DOUBLE_COMPLEX,pcol(2),mpi_comm_cols,mpierr)

   if(my_prow==prow(1) .and. my_pcol==pcol(1)) d(1) = a(1,1)

   deallocate(tmp, vr, ur, vc, uc, vur, uvc)

   ! distribute the arrays d and e to all processors

   allocate(tmpr(na))
   tmpr = d
   call mpi_allreduce(tmpr,d,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmpr = d
   call mpi_allreduce(tmpr,d,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   tmpr = e
   call mpi_allreduce(tmpr,e,na,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
   tmpr = e
   call mpi_allreduce(tmpr,e,na,MPI_REAL8,MPI_SUM,mpi_comm_cols,mpierr)
   deallocate(tmpr)

end subroutine tridiag_complex_2013

!-------------------------------------------------------------------------------

subroutine trans_ev_complex_2013(na, nqc, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  trans_ev_complex: Transforms the eigenvectors of a tridiagonal matrix back
!                    to the eigenvectors of the original matrix
!                    (like Scalapack Routine PZUNMTR)
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after tridiag_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tau(na)     Factors of the Householder vectors
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, nqc, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), q(ldq,*), tau(na)

   integer :: max_stored_rows

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, nstor
   integer istep, i, n, nc, ic, ics, ice, nb, cur_pcol

   complex*16, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)
   complex*16, allocatable:: tmat(:,:), h1(:), h2(:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   totalblocks = (na-1)/nblk + 1
   max_blocks_row = (totalblocks-1)/np_rows + 1
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk


   max_stored_rows = (63/nblk+1)*nblk

   allocate(tmat(max_stored_rows,max_stored_rows))
   allocate(h1(max_stored_rows*max_stored_rows))
   allocate(h2(max_stored_rows*max_stored_rows))
   allocate(tmp1(max_local_cols*max_stored_rows))
   allocate(tmp2(max_local_cols*max_stored_rows))
   allocate(hvb(max_local_rows*nblk))
   allocate(hvm(max_local_rows,max_stored_rows))

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

   l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

   nstor = 0

   ! In the complex case tau(2) /= 0
   if(my_prow == prow(1)) then
      q(1,1:l_cols) = q(1,1:l_cols)*((1.d0,0.d0)-tau(2))
   endif

   do istep=1,na,nblk

      ics = MAX(istep,3)
      ice = MIN(istep+nblk-1,na)
      if(ice<ics) cycle

      cur_pcol = pcol(istep)

      nb = 0
      do ic=ics,ice

         l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder vector
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector


         if(my_pcol==cur_pcol) then
            hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)
            if(my_prow==prow(ic-1)) then
               hvb(nb+l_rows) = 1.
            endif
         endif

         nb = nb+l_rows
      enddo

      if(nb>0) &
         call MPI_Bcast(hvb,nb,MPI_DOUBLE_COMPLEX,cur_pcol,mpi_comm_cols,mpierr)

      nb = 0
      do ic=ics,ice
         l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
         hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
         nstor = nstor+1
         nb = nb+l_rows
      enddo

      ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
      if(nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32)) then

         ! Calculate scalar products of stored vectors.
         ! This can be done in different ways, we use zherk

         tmat = 0
         if(l_rows>0) &
            call zherk('U','C',nstor,l_rows,CONE,hvm,ubound(hvm,1),CZERO,tmat,max_stored_rows)

         nc = 0
         do n=1,nstor-1
            h1(nc+1:nc+n) = tmat(1:n,n+1)
            nc = nc+n
         enddo

         if(nc>0) call mpi_allreduce(h1,h2,nc,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         ! Calculate triangular matrix T

         nc = 0
         tmat(1,1) = tau(ice-nstor+1)
         do n=1,nstor-1
            call ztrmv('L','C','N',n,tmat,max_stored_rows,h2(nc+1),1)
            tmat(n+1,1:n) = -conjg(h2(nc+1:nc+n))*tau(ice-nstor+n+1)
            tmat(n+1,n+1) = tau(ice-nstor+n+1)
            nc = nc+n
         enddo

         ! Q = Q - V * T * V**T * Q

         if(l_rows>0) then
            call zgemm('C','N',nstor,l_cols,l_rows,CONE,hvm,ubound(hvm,1), &
                       q,ldq,CZERO,tmp1,nstor)
         else
            tmp1(1:l_cols*nstor) = 0
         endif
         call mpi_allreduce(tmp1,tmp2,nstor*l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
         if(l_rows>0) then
            call ztrmm('L','L','N','N',nstor,l_cols,CONE,tmat,max_stored_rows,tmp2,nstor)
            call zgemm('N','N',l_rows,l_cols,nstor,-CONE,hvm,ubound(hvm,1), &
                       tmp2,nstor,CONE,q,ldq)
         endif
         nstor = 0
      endif

   enddo

   deallocate(tmat, h1, h2, tmp1, tmp2, hvb, hvm)


end subroutine trans_ev_complex_2013

!-------------------------------------------------------------------------------

subroutine mult_ah_b_complex_2013(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

!-------------------------------------------------------------------------------
!  mult_ah_b_complex:  Performs C := A**H * B
!
!      where:  A is a square matrix (na,na) which is optionally upper or lower triangular
!              B is a (na,ncb) matrix
!              C is a (na,ncb) matrix where optionally only the upper or lower
!              triangle may be computed
!
!  Parameters
!
!  uplo_a      'U' if A is upper triangular
!              'L' if A is lower triangular
!              anything else if A is a full matrix
!              Please note: This pertains to the original A (as set in the calling program)
!              whereas the transpose of A is used for calculations
!              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!              i.e. it may contain arbitrary numbers
!
!  uplo_c      'U' if only the upper diagonal part of C is needed
!              'L' if only the upper diagonal part of C is needed
!              anything else if the full matrix C is needed
!              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!              written to a certain extent, i.e. one shouldn't rely on the content there!
!
!  na          Number of rows/columns of A, number of rows of B and C
!
!  ncb         Number of columns  of B and C
!
!  a           Matrix A
!
!  lda         Leading dimension of a
!
!  b           Matrix B
!
!  ldb         Leading dimension of b
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  c           Matrix C
!
!  ldc         Leading dimension of c
!
!-------------------------------------------------------------------------------

   implicit none

   character*1 uplo_a, uplo_c

   integer na, ncb, lda, ldb, nblk, mpi_comm_rows, mpi_comm_cols, ldc
   complex*16 a(lda,*), b(ldb,*), c(ldc,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_rows_np
   integer np, n, nb, nblk_mult, lrs, lre, lcs, lce
   integer gcol_min, gcol, goff
   integer nstor, nr_done, noff, np_bc, n_aux_bc, nvals
   integer, allocatable :: lrs_save(:), lre_save(:)

   logical a_lower, a_upper, c_lower, c_upper

   complex*16, allocatable:: aux_mat(:,:), aux_bc(:), tmp1(:,:), tmp2(:,:)


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
   l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

   ! Block factor for matrix multiplications, must be a multiple of nblk

   if(na/np_rows<=256) then
      nblk_mult = (31/nblk+1)*nblk
   else
      nblk_mult = (63/nblk+1)*nblk
   endif

   allocate(aux_mat(l_rows,nblk_mult))
   allocate(aux_bc(l_rows*nblk))
   allocate(lrs_save(nblk))
   allocate(lre_save(nblk))

   a_lower = .false.
   a_upper = .false.
   c_lower = .false.
   c_upper = .false.

   if(uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
   if(uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
   if(uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
   if(uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

   ! Build up the result matrix by processor rows

   do np = 0, np_rows-1

      ! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0   ! Number of columns stored in aux_mat

      ! Loop over the blocks on row np

      do nb=0,(l_rows_np-1)/nblk

         goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

         ! Get the processor column which owns this block (A is transposed, so we need the column)
         ! and the offset in blocks within this column.
         ! The corresponding block column in A is then broadcast to all for multiplication with B

         np_bc = MOD(goff,np_cols)
         noff = goff/np_cols
         n_aux_bc = 0

         ! Gather up the complete block column of A on the owner

         do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

            gcol = goff*nblk + n ! global column corresponding to n
            if(nstor==0 .and. n==1) gcol_min = gcol

            lrs = 1       ! 1st local row number for broadcast
            lre = l_rows  ! last local row number for broadcast
            if(a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
            if(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            if(lrs<=lre) then
               nvals = lre-lrs+1
               if(my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
               n_aux_bc = n_aux_bc + nvals
            endif

            lrs_save(n) = lrs
            lre_save(n) = lre

         enddo

         ! Broadcast block column

         call MPI_Bcast(aux_bc,n_aux_bc,MPI_DOUBLE_COMPLEX,np_bc,mpi_comm_cols,mpierr)

         ! Insert what we got in aux_mat

         n_aux_bc = 0
         do n = 1, min(l_rows_np-nb*nblk,nblk)
            nstor = nstor+1
            lrs = lrs_save(n)
            lre = lre_save(n)
            if(lrs<=lre) then
               nvals = lre-lrs+1
               aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
               n_aux_bc = n_aux_bc + nvals
            endif
         enddo

         ! If we got nblk_mult columns in aux_mat or this is the last block
         ! do the matrix multiplication

         if(nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

            lrs = 1       ! 1st local row number for multiply
            lre = l_rows  ! last local row number for multiply
            if(a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
            if(a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

            lcs = 1       ! 1st local col number for multiply
            lce = l_cols  ! last local col number for multiply
            if(c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
            if(c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

            if(lcs<=lce) then
               allocate(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce))
               if(lrs<=lre) then
                  call zgemm('C','N',nstor,lce-lcs+1,lre-lrs+1,(1.d0,0.d0),aux_mat(lrs,1),ubound(aux_mat,1), &
                             b(lrs,lcs),ldb,(0.d0,0.d0),tmp1,nstor)
               else
                  tmp1 = 0
               endif

               ! Sum up the results and send to processor row np
               call mpi_reduce(tmp1,tmp2,nstor*(lce-lcs+1),MPI_DOUBLE_COMPLEX,MPI_SUM,np,mpi_comm_rows,mpierr)

               ! Put the result into C
               if(my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

               deallocate(tmp1,tmp2)
            endif

            nr_done = nr_done+nstor
            nstor=0
            aux_mat(:,:)=0
         endif
      enddo
   enddo

   deallocate(aux_mat, aux_bc, lrs_save, lre_save)

end subroutine mult_ah_b_complex_2013

!-------------------------------------------------------------------------------

subroutine solve_tridi_2013( na, nev, d, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols )

   implicit none

   integer  na, nev, ldq, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 d(na), e(na), q(ldq,*)

   integer i, j, n, np, nc, nev1, l_cols, l_rows
   integer my_prow, my_pcol, np_rows, np_cols, mpierr

   integer, allocatable :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

   ! Set Q to 0

   q(1:l_rows, 1:l_cols) = 0.

   ! Get the limits of the subdivisons, each subdivison has as many cols
   ! as fit on the respective processor column

   allocate(limits(0:np_cols))

   limits(0) = 0
   do np=0,np_cols-1
      nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

      ! Check for the case that a column has have zero width.
      ! This is not supported!
      ! Scalapack supports it but delivers no results for these columns,
      ! which is rather annoying
      if(nc==0) then
         print *,'ERROR: Problem contains processor column with zero width'
         call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
      endif

      limits(np+1) = limits(np) + nc
   enddo

   ! Subdivide matrix by subtracting rank 1 modifications

   do i=1,np_cols-1
      n = limits(i)
      d(n) = d(n)-abs(e(n))
      d(n+1) = d(n+1)-abs(e(n))
   enddo

   ! Solve sub problems on processsor columns

   nc = limits(my_pcol) ! column after which my problem starts

   if(np_cols>1) then
      nev1 = l_cols ! all eigenvectors are needed
   else
      nev1 = MIN(nev,l_cols)
   endif
   call solve_tridi_col(l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk, mpi_comm_rows)


   ! If there is only 1 processor column, we are done

   if(np_cols==1) then
      deallocate(limits)
      return
   endif

   ! Set index arrays for Q columns

   ! Dense distribution scheme:

   allocate(l_col(na))
   allocate(p_col(na))

   n = 0
   do np=0,np_cols-1
      nc = local_index(na, np, np_cols, nblk, -1)
      do i=1,nc
         n = n+1
         l_col(n) = i
         p_col(n) = np
      enddo
   enddo

   ! Block cyclic distribution scheme, only nev columns are set:

   allocate(l_col_bc(na))
   allocate(p_col_bc(na))
   p_col_bc(:) = -1
   l_col_bc(:) = -1

   do i = 0, na-1, nblk*np_cols
      do j = 0, np_cols-1
         do n = 1, nblk
            if(i+j*nblk+n <= MIN(nev,na)) then
               p_col_bc(i+j*nblk+n) = j
               l_col_bc(i+j*nblk+n) = i/np_cols + n
            endif
         enddo
      enddo
   enddo

   ! Recursively merge sub problems

   call merge_recursive(0, np_cols)

   deallocate(limits,l_col,p_col,l_col_bc,p_col_bc)

contains
recursive subroutine merge_recursive(np_off, nprocs)

   implicit none

   ! noff is always a multiple of nblk_ev
   ! nlen-noff is always > nblk_ev

   integer np_off, nprocs
   integer np1, np2, noff, nlen, nmid, n, mpi_status(mpi_status_size)

   if(nprocs<=1) then
      ! Safety check only
      print *,"INTERNAL error merge_recursive: nprocs=",nprocs
      call mpi_abort(MPI_COMM_WORLD,1,mpierr)
   endif

   ! Split problem into 2 subproblems of size np1 / np2

   np1 = nprocs/2
   np2 = nprocs-np1

   if(np1 > 1) call merge_recursive(np_off, np1)
   if(np2 > 1) call merge_recursive(np_off+np1, np2)


   noff = limits(np_off)
   nmid = limits(np_off+np1) - noff
   nlen = limits(np_off+nprocs) - noff

   if(my_pcol==np_off) then
      do n=np_off+np1,np_off+nprocs-1
         call mpi_send(d(noff+1),nmid,MPI_REAL8,n,1,mpi_comm_cols,mpierr)
      enddo
   endif
   if(my_pcol>=np_off+np1 .and. my_pcol<np_off+nprocs) then
      call mpi_recv(d(noff+1),nmid,MPI_REAL8,np_off,1,mpi_comm_cols,mpi_status,mpierr)
   endif

   if(my_pcol==np_off+np1) then
      do n=np_off,np_off+np1-1
         call mpi_send(d(noff+nmid+1),nlen-nmid,MPI_REAL8,n,1,mpi_comm_cols,mpierr)
      enddo
   endif
   if(my_pcol>=np_off .and. my_pcol<np_off+np1) then
      call mpi_recv(d(noff+nmid+1),nlen-nmid,MPI_REAL8,np_off+np1,1,mpi_comm_cols,mpi_status,mpierr)
   endif

   if(nprocs == np_cols) then

      ! Last merge, result distribution must be block cyclic, noff==0,
      ! p_col_bc is set so that only nev eigenvalues are calculated

      call merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
                         nblk, mpi_comm_rows, mpi_comm_cols, l_col, p_col, &
                         l_col_bc, p_col_bc, np_off, nprocs )

   else

      ! Not last merge, leave dense column distribution

      call merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
                         nblk, mpi_comm_rows, mpi_comm_cols, l_col(noff+1), p_col(noff+1), &
                         l_col(noff+1), p_col(noff+1), np_off, nprocs )
   endif

end subroutine merge_recursive

end subroutine solve_tridi_2013

!-------------------------------------------------------------------------------

subroutine solve_tridi_col( na, nev, nqoff, d, e, q, ldq, nblk, mpi_comm_rows )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!

   implicit none

   integer  na, nev, nqoff, ldq, nblk, mpi_comm_rows
   real*8 d(na), e(na), q(ldq,*)

   integer, parameter:: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

   real*8, allocatable :: qmat1(:,:), qmat2(:,:)

   integer i, n, np
   integer ndiv, noff, nmid, nlen, max_size
   integer my_prow, np_rows, mpierr

   integer, allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)

   ! Calculate the number of subdivisions needed.

   n = na
   ndiv = 1
   do while(2*ndiv<=np_rows .and. n>2*min_submatrix_size)
      n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
      ndiv = ndiv*2
   enddo

   ! If there is only 1 processor row and not all eigenvectors are needed
   ! and the matrix size is big enough, then use 2 subdivisions
   ! so that merge_systems is called once and only the needed
   ! eigenvectors are calculated for the final problem.

   if(np_rows==1 .and. nev<na .and. na>2*min_submatrix_size) ndiv = 2

   allocate(limits(0:ndiv))

   limits(0) = 0
   limits(ndiv) = na

   n = ndiv
   do while(n>1)
      n = n/2 ! n is always a power of 2
      do i=0,ndiv-1,2*n
         ! We want to have even boundaries (for cache line alignments)
         limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
      enddo
   enddo

   ! Calculate the maximum size of a subproblem

   max_size = 0
   do i=1,ndiv
      max_size = MAX(max_size,limits(i)-limits(i-1))
   enddo

   ! Subdivide matrix by subtracting rank 1 modifications

   do i=1,ndiv-1
      n = limits(i)
      d(n) = d(n)-abs(e(n))
      d(n+1) = d(n+1)-abs(e(n))
   enddo

   if(np_rows==1)    then

      ! For 1 processor row there may be 1 or 2 subdivisions

      do n=0,ndiv-1
         noff = limits(n)        ! Start of subproblem
         nlen = limits(n+1)-noff ! Size of subproblem

         call solve_tridi_single(nlen,d(noff+1),e(noff+1),q(nqoff+noff+1,noff+1),ubound(q,1))
      enddo

   else

      ! Solve sub problems in parallel with solve_tridi_single
      ! There is at maximum 1 subproblem per processor

      allocate(qmat1(max_size,max_size))
      allocate(qmat2(max_size,max_size))

      qmat1 = 0 ! Make sure that all elements are defined

      if(my_prow < ndiv) then

         noff = limits(my_prow)        ! Start of subproblem
         nlen = limits(my_prow+1)-noff ! Size of subproblem

         call solve_tridi_single(nlen,d(noff+1),e(noff+1),qmat1,ubound(qmat1,1))

      endif

      ! Fill eigenvectors in qmat1 into global matrix q

      do np = 0, ndiv-1

         noff = limits(np)
         nlen = limits(np+1)-noff

         call MPI_Bcast(d(noff+1),nlen,MPI_REAL8,np,mpi_comm_rows,mpierr)
         qmat2 = qmat1
         call MPI_Bcast(qmat2,max_size*max_size,MPI_REAL8,np,mpi_comm_rows,mpierr)

         do i=1,nlen
            call distribute_global_column(qmat2(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
         enddo

      enddo

      deallocate(qmat1, qmat2)

   endif

   ! Allocate and set index arrays l_col and p_col

   allocate(l_col(na), p_col_i(na),  p_col_o(na))

   do i=1,na
      l_col(i) = i
      p_col_i(i) = 0
      p_col_o(i) = 0
   enddo

   ! Merge subproblems

   n = 1
   do while(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

      do i=0,ndiv-1,2*n

         noff = limits(i)
         nmid = limits(i+n) - noff
         nlen = limits(i+2*n) - noff

        if(nlen == na) then
           ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
           p_col_o(nev+1:na) = -1
        endif

         call merge_systems(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                            mpi_comm_rows, mpi_comm_self, l_col(noff+1), p_col_i(noff+1), &
                            l_col(noff+1), p_col_o(noff+1), 0, 1)

      enddo

      n = 2*n

   enddo

   deallocate(limits, l_col, p_col_i, p_col_o)

end subroutine solve_tridi_col

!-------------------------------------------------------------------------------

subroutine solve_tridi_single(nlen, d, e, q, ldq)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.

   implicit none

   integer nlen, ldq
   real*8 d(nlen), e(nlen), q(ldq,nlen)

   real*8, allocatable :: work(:), qtmp(:), ds(:), es(:)
   real*8 dtmp

   integer i, j, lwork, liwork, info, mpierr
   integer, allocatable :: iwork(:)

   allocate(ds(nlen), es(nlen))

   ! Save d and e for the case that dstedc fails

   ds(:) = d(:)
   es(:) = e(:)

   ! First try dstedc, this is normally faster but it may fail sometimes (why???)

   lwork = 1 + 4*nlen + nlen**2
   liwork =  3 + 5*nlen
   allocate(work(lwork), iwork(liwork))
   call dstedc('I',nlen,d,e,q,ldq,work,lwork,iwork,liwork,info)

   if(info /= 0) then

      ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

      print '(a,i8,a)','Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

      d(:) = ds(:)
      e(:) = es(:)
      call dsteqr('I',nlen,d,e,q,ldq,work,info)

      ! If DSTEQR fails also, we don't know what to do further ...
      if(info /= 0) then
         print '(a,i8,a)','ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif

   end if

   deallocate(work,iwork,ds,es)

   ! Check if eigenvalues are monotonically increasing
   ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

   do i=1,nlen-1
      if(d(i+1)<d(i)) then
         if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1d-14) then
            print '(a,i8,2g25.16)','***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
! XR: not really useful, and become annoying for SOSEX calculations; comment out for now
!         else
!            print '(a,i8,2g25.16)','Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
!            print '(a)', 'The eigenvalues from a lapack call are not sorted to machine precision.'
!            print '(a)', 'In this extent, this is completely harmless.'
!            print '(a)', 'Still, we keep this info message just in case.'
         end if
         allocate(qtmp(nlen))
         dtmp = d(i+1)
         qtmp(1:nlen) = q(1:nlen,i+1)
         do j=i,1,-1
            if(dtmp<d(j)) then
               d(j+1)        = d(j)
               q(1:nlen,j+1) = q(1:nlen,j)
            else
               exit ! Loop
            endif
         enddo
         d(j+1)        = dtmp
         q(1:nlen,j+1) = qtmp(1:nlen)
         deallocate(qtmp)
      endif
   enddo

end subroutine solve_tridi_single

!-------------------------------------------------------------------------------

subroutine merge_systems( na, nm, d, e, q, ldq, nqoff, nblk, mpi_comm_rows, mpi_comm_cols, &
                          l_col, p_col, l_col_out, p_col_out, npc_0, npc_n)

   implicit none

   integer  na, nm, ldq, nqoff, nblk, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n
   integer  l_col(na), p_col(na), l_col_out(na), p_col_out(na)
   real*8 d(na), e, q(ldq,*)

   integer, parameter :: max_strip=128

   real*8 beta, sig, s, c, t, tau, rho, eps, tol, dlamch, dlapy2, qtrans(2,2), dmax, zmax, d1new, d2new
   real*8 z(na), d1(na), d2(na), z1(na), delta(na), dbase(na), ddiff(na), ev_scale(na), tmp(na)
   real*8 d1u(na), zu(na), d1l(na), zl(na)
   real*8, allocatable :: qtmp1(:,:), qtmp2(:,:), ev(:,:)

   real*8, allocatable :: z_p(:,:)


   integer i, j, na1, na2, l_rows, l_cols, l_rqs, l_rqe, l_rqm, ns, info
   integer l_rnm, nnzu, nnzl, ndef, ncnt, max_local_cols, l_cols_qreorg, np, l_idx, nqcols1, nqcols2
   integer my_proc, n_procs, my_prow, my_pcol, np_rows, np_cols, mpierr, mpi_status(mpi_status_size)
   integer np_next, np_prev, np_rem
   integer idx(na), idx1(na), idx2(na)
   integer coltyp(na), idxq1(na), idxq2(na)

   integer max_threads, my_thread
   integer omp_get_max_threads, omp_get_thread_num

   max_threads = omp_get_max_threads()

   allocate(z_p(na,0:max_threads-1))

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! If my processor column isn't in the requested set, do nothing

   if(my_pcol<npc_0 .or. my_pcol>=npc_0+npc_n) return

   ! Determine number of "next" and "prev" column for ring sends

   if(my_pcol == npc_0+npc_n-1) then
      np_next = npc_0
   else
      np_next = my_pcol + 1
   endif

   if(my_pcol == npc_0) then
      np_prev = npc_0+npc_n-1
   else
      np_prev = my_pcol - 1
   endif

   call check_monotony(nm,d,'Input1')
   call check_monotony(na-nm,d(nm+1),'Input2')

   ! Get global number of processors and my processor number.
   ! Please note that my_proc does not need to match any real processor number,
   ! it is just used for load balancing some loops.

   n_procs = np_rows*npc_n
   my_proc = my_prow*npc_n + (my_pcol-npc_0) ! Row major


   ! Local limits of the rows of Q

   l_rqs = local_index(nqoff+1 , my_prow, np_rows, nblk, +1) ! First row of Q
   l_rqm = local_index(nqoff+nm, my_prow, np_rows, nblk, -1) ! Last row <= nm
   l_rqe = local_index(nqoff+na, my_prow, np_rows, nblk, -1) ! Last row of Q

   l_rnm  = l_rqm-l_rqs+1 ! Number of local rows <= nm
   l_rows = l_rqe-l_rqs+1 ! Total number of local rows


   ! My number of local columns

   l_cols = COUNT(p_col(1:na)==my_pcol)

   ! Get max number of local columns

   max_local_cols = 0
   do np = npc_0, npc_0+npc_n-1
      max_local_cols = MAX(max_local_cols,COUNT(p_col(1:na)==np))
   enddo



   ! Calculations start here

   beta = abs(e)
   sig  = sign(1.d0,e)

   ! Calculate rank-1 modifier z

   z(:) = 0

   if(MOD((nqoff+nm-1)/nblk,np_rows)==my_prow) then
      ! nm is local on my row
      do i = 1, na
         if(p_col(i)==my_pcol) z(i) = q(l_rqm,l_col(i))
      enddo
   endif

   if(MOD((nqoff+nm)/nblk,np_rows)==my_prow) then
      ! nm+1 is local on my row
      do i = 1, na
         if(p_col(i)==my_pcol) z(i) = z(i) + sig*q(l_rqm+1,l_col(i))
      enddo
   endif

   call global_gather(z, na)

   ! Normalize z so that norm(z) = 1.  Since z is the concatenation of
   ! two normalized vectors, norm2(z) = sqrt(2).

   z = z/sqrt(2.0d0)
   rho = 2.*beta

   ! Calculate index for merging both systems by ascending eigenvalues

   call DLAMRG( nm, na-nm, d, 1, 1, idx )

   ! Calculate the allowable deflation tolerance

   zmax = maxval(abs(z))
   dmax = maxval(abs(d))
   EPS = DLAMCH( 'Epsilon' )
   TOL = 8.*EPS*MAX(dmax,zmax)

   ! If the rank-1 modifier is small enough, no more needs to be done
   ! except to reorganize D and Q

   IF( RHO*zmax <= TOL ) THEN

      ! Rearrange eigenvalues

      tmp = d
      do i=1,na
         d(i) = tmp(idx(i))
      enddo

      ! Rearrange eigenvectors

      call resort_ev(idx)

      return
   ENDIF

   ! Merge and deflate system

   na1 = 0
   na2 = 0

   ! COLTYP:
   ! 1 : non-zero in the upper half only;
   ! 2 : dense;
   ! 3 : non-zero in the lower half only;
   ! 4 : deflated.

   coltyp(1:nm) = 1
   coltyp(nm+1:na) = 3

   do i=1,na

      if(rho*abs(z(idx(i))) <= tol) then

         ! Deflate due to small z component.

         na2 = na2+1
         d2(na2)   = d(idx(i))
         idx2(na2) = idx(i)
         coltyp(idx(i)) = 4

      else if(na1>0) then

         ! Check if eigenvalues are close enough to allow deflation.

         S = Z(idx(i))
         C = Z1(na1)

         ! Find sqrt(a**2+b**2) without overflow or
         ! destructive underflow.

         TAU = DLAPY2( C, S )
         T = D1(na1) - D(idx(i))
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ) <= TOL ) THEN

            ! Deflation is possible.

            na2 = na2+1

            Z1(na1) = TAU

            d2new = D(idx(i))*C**2 + D1(na1)*S**2
            d1new = D(idx(i))*S**2 + D1(na1)*C**2

            ! D(idx(i)) >= D1(na1) and C**2 + S**2 == 1.0
            ! This means that after the above transformation it must be
            !    D1(na1) <= d1new <= D(idx(i))
            !    D1(na1) <= d2new <= D(idx(i))
            !
            ! D1(na1) may get bigger but it is still smaller than the next D(idx(i+1))
            ! so there is no problem with sorting here.
            ! d2new <= D(idx(i)) which means that it might be smaller than D2(na2-1)
            ! which makes a check (and possibly a resort) necessary.
            !
            ! The above relations may not hold exactly due to numeric differences
            ! so they have to be enforced in order not to get troubles with sorting.


            if(d1new<D1(na1)  ) d1new = D1(na1)
            if(d1new>D(idx(i))) d1new = D(idx(i))

            if(d2new<D1(na1)  ) d2new = D1(na1)
            if(d2new>D(idx(i))) d2new = D(idx(i))

            D1(na1) = d1new

            do j=na2-1,1,-1
               if(d2new<d2(j)) then
                  d2(j+1)   = d2(j)
                  idx2(j+1) = idx2(j)
               else
                  exit ! Loop
               endif
            enddo

            d2(j+1)   = d2new
            idx2(j+1) = idx(i)

            qtrans(1,1) = C; qtrans(1,2) =-S
            qtrans(2,1) = S; qtrans(2,2) = C

            call transform_columns(idx(i), idx1(na1))

            if(coltyp(idx(i))==1 .and. coltyp(idx1(na1))/=1) coltyp(idx1(na1)) = 2
            if(coltyp(idx(i))==3 .and. coltyp(idx1(na1))/=3) coltyp(idx1(na1)) = 2

            coltyp(idx(i)) = 4

         else
            na1 = na1+1
            d1(na1) = d(idx(i))
            z1(na1) = z(idx(i))
            idx1(na1) = idx(i)
         endif
      else
         na1 = na1+1
         d1(na1) = d(idx(i))
         z1(na1) = z(idx(i))
         idx1(na1) = idx(i)
      endif

   enddo
   call check_monotony(na1,d1,'Sorted1')
   call check_monotony(na2,d2,'Sorted2')

   if(na1==1 .or. na1==2) then
      ! if(my_proc==0) print *,'--- Remark solve_tridi: na1==',na1,' proc==',myid

      if(na1==1) then
         d(1) = d1(1) + rho*z1(1)**2 ! solve secular equation
      else ! na1==2
         call DLAED5(1, d1, z1, qtrans(1,1), rho, d(1))
         call DLAED5(2, d1, z1, qtrans(1,2), rho, d(2))

         call transform_columns(idx1(1), idx1(2))
      endif

      ! Add the deflated eigenvalues
      d(na1+1:na) = d2(1:na2)

      ! Calculate arrangement of all eigenvalues  in output

      call DLAMRG( na1, na-na1, d, 1, 1, idx )

      ! Rearrange eigenvalues

      tmp = d
      do i=1,na
         d(i) = tmp(idx(i))
      enddo

      ! Rearrange eigenvectors

      do i=1,na
         if(idx(i)<=na1) then
            idxq1(i) = idx1(idx(i))
         else
            idxq1(i) = idx2(idx(i)-na1)
         endif
      enddo

      call resort_ev(idxq1)

   else if(na1>2) then

      ! Solve secular equation

      z(1:na1) = 1

      z_p(1:na1,:) = 1

      dbase(1:na1) = 0
      ddiff(1:na1) = 0

      info = 0

!$OMP PARALLEL PRIVATE(i,my_thread,delta,s,info,j)
      my_thread = omp_get_thread_num()
!$OMP DO

      DO i = my_proc+1, na1, n_procs ! work distributed over all processors

         call DLAED4_OWN(na1, i, d1, z1, delta, rho, s, info) ! s is not used!

         if(info/=0) then
            ! If DLAED4 fails (may happen especially for LAPACK versions before 3.2)
            ! use the more stable bisection algorithm in solve_secular_equation
            ! print *,'ERROR DLAED4 n=',na1,'i=',i,' Using Bisection'
            call solve_secular_equation(na1, i, d1, z1, delta, rho, s)
         endif

         ! Compute updated z

         do j=1,na1
            if(i/=j)  z_p(j,my_thread) = z_p(j,my_thread)*( delta(j) / (d1(j)-d1(i)) )
         enddo
         z_p(i,my_thread) = z_p(i,my_thread)*delta(i)

         ! store dbase/ddiff

         if(i<na1) then
            if(abs(delta(i+1)) < abs(delta(i))) then
               dbase(i) = d1(i+1)
               ddiff(i) = delta(i+1)
            else
               dbase(i) = d1(i)
               ddiff(i) = delta(i)
            endif
         else
            dbase(i) = d1(i)
            ddiff(i) = delta(i)
         endif
      enddo

!$OMP END PARALLEL

      do i = 0, max_threads-1
         z(1:na1) = z(1:na1)*z_p(1:na1,i)
      enddo


      call global_product(z, na1)
      z(1:na1) = SIGN( SQRT( -z(1:na1) ), z1(1:na1) )

      call global_gather(dbase, na1)
      call global_gather(ddiff, na1)
      d(1:na1) = dbase(1:na1) - ddiff(1:na1)

      ! Calculate scale factors for eigenvectors

      ev_scale(:) = 0.

!$OMP PARALLEL DO PRIVATE(i) SHARED(na1, my_proc, n_procs,  &
!$OMP d1,dbase, ddiff, z, ev_scale) &
!$OMP DEFAULT(NONE)

      DO i = my_proc+1, na1, n_procs ! work distributed over all processors

         ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
         ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

         ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
         ! in exactly this order, but we want to prevent compiler optimization
!         ev_scale_val = ev_scale(i)
         call add_tmp(d1, dbase, ddiff, z, ev_scale(i), na1,i)
!         ev_scale(i) = ev_scale_val
      enddo

!$OMP END PARALLEL DO


      call global_gather(ev_scale, na1)

      ! Add the deflated eigenvalues
      d(na1+1:na) = d2(1:na2)

      ! Calculate arrangement of all eigenvalues  in output

      call DLAMRG( na1, na-na1, d, 1, 1, idx )

      ! Rearrange eigenvalues

      tmp = d
      do i=1,na
         d(i) = tmp(idx(i))
      enddo
      call check_monotony(na,d,'Output')

      ! Eigenvector calculations


      ! Calculate the number of columns in the new local matrix Q
      ! which are updated from non-deflated/deflated eigenvectors.
      ! idxq1/2 stores the global column numbers.

      nqcols1 = 0 ! number of non-deflated eigenvectors
      nqcols2 = 0 ! number of deflated eigenvectors
      DO i = 1, na
         if(p_col_out(i)==my_pcol) then
            if(idx(i)<=na1) then
               nqcols1 = nqcols1+1
               idxq1(nqcols1) = i
            else
               nqcols2 = nqcols2+1
               idxq2(nqcols2) = i
            endif
         endif
      enddo

      allocate(ev(max_local_cols,MIN(max_strip,MAX(1,nqcols1))))
      allocate(qtmp1(MAX(1,l_rows),max_local_cols))
      allocate(qtmp2(MAX(1,l_rows),MIN(max_strip,MAX(1,nqcols1))))

      ! Gather nonzero upper/lower components of old matrix Q
      ! which are needed for multiplication with new eigenvectors

      qtmp1 = 0 ! May contain empty (unset) parts
      qtmp2 = 0 ! Not really needed

      nnzu = 0
      nnzl = 0
      do i = 1, na1
         l_idx = l_col(idx1(i))
         if(p_col(idx1(i))==my_pcol) then
            if(coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
               nnzu = nnzu+1
               qtmp1(1:l_rnm,nnzu) = q(l_rqs:l_rqm,l_idx)
            endif
            if(coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
               nnzl = nnzl+1
               qtmp1(l_rnm+1:l_rows,nnzl) = q(l_rqm+1:l_rqe,l_idx)
            endif
         endif
      enddo

      ! Gather deflated eigenvalues behind nonzero components

      ndef = max(nnzu,nnzl)
      do i = 1, na2
         l_idx = l_col(idx2(i))
         if(p_col(idx2(i))==my_pcol) then
            ndef = ndef+1
            qtmp1(1:l_rows,ndef) = q(l_rqs:l_rqe,l_idx)
         endif
      enddo

      l_cols_qreorg = ndef ! Number of columns in reorganized matrix

      ! Set (output) Q to 0, it will sum up new Q

      DO i = 1, na
         if(p_col_out(i)==my_pcol) q(l_rqs:l_rqe,l_col_out(i)) = 0
      enddo


      np_rem = my_pcol

      do np = 1, npc_n

         ! Do a ring send of qtmp1

         if(np>1) then

            if(np_rem==npc_0) then
               np_rem = npc_0+npc_n-1
            else
               np_rem = np_rem-1
            endif

            call MPI_Sendrecv_replace(qtmp1, l_rows*max_local_cols, MPI_REAL8, &
                                      np_next, 1111, np_prev, 1111, &
                                      mpi_comm_cols, mpi_status, mpierr)
         endif

         ! Gather the parts in d1 and z which are fitting to qtmp1.
         ! This also delivers nnzu/nnzl for proc np_rem

         nnzu = 0
         nnzl = 0
         do i=1,na1
            if(p_col(idx1(i))==np_rem) then
               if(coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
                  nnzu = nnzu+1
                  d1u(nnzu) = d1(i)
                  zu (nnzu) = z (i)
               endif
               if(coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
                  nnzl = nnzl+1
                  d1l(nnzl) = d1(i)
                  zl (nnzl) = z (i)
               endif
            endif
         enddo

         ! Set the deflated eigenvectors in Q (comming from proc np_rem)

         ndef = MAX(nnzu,nnzl) ! Remote counter in input matrix
         do i = 1, na
            j = idx(i)
            if(j>na1) then
               if(p_col(idx2(j-na1))==np_rem) then
                  ndef = ndef+1
                  if(p_col_out(i)==my_pcol) &
                     q(l_rqs:l_rqe,l_col_out(i)) = qtmp1(1:l_rows,ndef)
               endif
            endif
         enddo

         do ns = 0, nqcols1-1, max_strip ! strimining loop

            ncnt = MIN(max_strip,nqcols1-ns) ! number of columns in this strip

            ! Get partial result from (output) Q

            do i = 1, ncnt
               qtmp2(1:l_rows,i) = q(l_rqs:l_rqe,l_col_out(idxq1(i+ns)))
            enddo

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with upper half of Q:

            do i = 1, ncnt
               j = idx(idxq1(i+ns))
               ! Calculate the j-th eigenvector of the deflated system
               ! See above why we are doing it this way!
               tmp(1:nnzu) = d1u(1:nnzu)-dbase(j)
               call v_add_s(tmp,nnzu,ddiff(j))
               ev(1:nnzu,i) = zu(1:nnzu) / tmp(1:nnzu) * ev_scale(j)
            enddo

            ! Multiply old Q with eigenvectors (upper half)

            if(l_rnm>0 .and. ncnt>0 .and. nnzu>0) &
               call dgemm('N','N',l_rnm,ncnt,nnzu,1.d0,qtmp1,ubound(qtmp1,1),ev,ubound(ev,1), &
                          1.d0,qtmp2(1,1),ubound(qtmp2,1))

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with lower half of Q:

            do i = 1, ncnt
               j = idx(idxq1(i+ns))
               ! Calculate the j-th eigenvector of the deflated system
               ! See above why we are doing it this way!
               tmp(1:nnzl) = d1l(1:nnzl)-dbase(j)
               call v_add_s(tmp,nnzl,ddiff(j))
               ev(1:nnzl,i) = zl(1:nnzl) / tmp(1:nnzl) * ev_scale(j)
            enddo

            ! Multiply old Q with eigenvectors (lower half)

            if(l_rows-l_rnm>0 .and. ncnt>0 .and. nnzl>0) &
               call dgemm('N','N',l_rows-l_rnm,ncnt,nnzl,1.d0,qtmp1(l_rnm+1,1),ubound(qtmp1,1),ev,ubound(ev,1), &
                          1.d0,qtmp2(l_rnm+1,1),ubound(qtmp2,1))

            ! Put partial result into (output) Q

            do i = 1, ncnt
               q(l_rqs:l_rqe,l_col_out(idxq1(i+ns))) = qtmp2(1:l_rows,i)
            enddo

         enddo
      enddo

      deallocate(ev, qtmp1, qtmp2)

   endif

!-------------------------------------------------------------------------------

contains
  subroutine add_tmp(d1, dbase, ddiff, z, ev_scale_value, na1,i)

    implicit none
    
    integer, intent(in) :: na1, i
    
    real*8, intent(in)  :: d1(:), dbase(:), ddiff(:), z(:)
    real*8, intent(inout) :: ev_scale_value
    real*8  :: tmp(1:na1)

         ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
         ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

         ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
         ! in exactly this order, but we want to prevent compiler optimization
    
    tmp(1:na1) = d1(1:na1) -dbase(i)
    call v_add_s(tmp(1:na1),na1,ddiff(i))
    
    tmp(1:na1) = z(1:na1) / tmp(1:na1)
   
    ev_scale_value = 1.0/sqrt(dot_product(tmp(1:na1),tmp(1:na1)))
    
  end subroutine add_tmp

subroutine resort_ev(idx_ev)

   implicit none

   integer idx_ev(*)
   integer i, nc, pc1, pc2, lc1, lc2, l_cols_out

   real*8, allocatable :: qtmp(:,:)

   if(l_rows==0) return ! My processor column has no work to do

   ! Resorts eigenvectors so that q_new(:,i) = q_old(:,idx_ev(i))

   l_cols_out = COUNT(p_col_out(1:na)==my_pcol)
   allocate(qtmp(l_rows,l_cols_out))


   nc = 0

   do i=1,na

      pc1 = p_col(idx_ev(i))
      lc1 = l_col(idx_ev(i))
      pc2 = p_col_out(i)

      if(pc2<0) cycle ! This column is not needed in output

      if(pc2==my_pcol) nc = nc+1 ! Counter for output columns

      if(pc1==my_pcol) then
         if(pc2==my_pcol) then
            ! send and recieve column are local
            qtmp(1:l_rows,nc) = q(l_rqs:l_rqe,lc1)
         else
            call mpi_send(q(l_rqs,lc1),l_rows,MPI_REAL8,pc2,mod(i,4096),mpi_comm_cols,mpierr)
         endif
      else if(pc2==my_pcol) then
         call mpi_recv(qtmp(1,nc),l_rows,MPI_REAL8,pc1,mod(i,4096),mpi_comm_cols,mpi_status,mpierr)
      endif
   enddo

   ! Insert qtmp into (output) q

   nc = 0

   do i=1,na

      pc2 = p_col_out(i)
      lc2 = l_col_out(i)

      if(pc2==my_pcol) then
         nc = nc+1
         q(l_rqs:l_rqe,lc2) = qtmp(1:l_rows,nc)
      endif
   enddo

   deallocate(qtmp)

end subroutine resort_ev

subroutine transform_columns(col1, col2)

   implicit none

   integer col1, col2
   integer pc1, pc2, lc1, lc2

   if(l_rows==0) return ! My processor column has no work to do

   pc1 = p_col(col1)
   lc1 = l_col(col1)
   pc2 = p_col(col2)
   lc2 = l_col(col2)

   if(pc1==my_pcol) then
      if(pc2==my_pcol) then
         ! both columns are local
         tmp(1:l_rows)      = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + q(l_rqs:l_rqe,lc2)*qtrans(2,1)
         q(l_rqs:l_rqe,lc2) = q(l_rqs:l_rqe,lc1)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
         q(l_rqs:l_rqe,lc1) = tmp(1:l_rows)
      else
         call mpi_sendrecv(q(l_rqs,lc1),l_rows,MPI_REAL8,pc2,1, &
                           tmp,l_rows,MPI_REAL8,pc2,1, &
                           mpi_comm_cols,mpi_status,mpierr)
         q(l_rqs:l_rqe,lc1) = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + tmp(1:l_rows)*qtrans(2,1)
      endif
   else if(pc2==my_pcol) then
      call mpi_sendrecv(q(l_rqs,lc2),l_rows,MPI_REAL8,pc1,1, &
                        tmp,l_rows,MPI_REAL8,pc1,1, &
                        mpi_comm_cols,mpi_status,mpierr)
      q(l_rqs:l_rqe,lc2) = tmp(1:l_rows)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
   endif

end subroutine transform_columns

subroutine global_gather(z, n)

   ! This routine sums up z over all processors.
   ! It should only be used for gathering distributed results,
   ! i.e. z(i) should be nonzero on exactly 1 processor column,
   ! otherways the results may be numerically different on different columns

   implicit none

   integer n
   real*8 z(n)

   real*8 tmp(n)

   if(npc_n==1 .and. np_rows==1) return ! nothing to do

   ! Do an mpi_allreduce over processor rows

   call mpi_allreduce(z, tmp, n, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

   ! If only 1 processor column, we are done
   if(npc_n==1) then
      z(:) = tmp(:)
      return
   endif

   ! If all processor columns are involved, we can use mpi_allreduce
   if(npc_n==np_cols) then
      call mpi_allreduce(tmp, z, n, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)
      return
   endif

   ! Do a ring send over processor columns
   z(:) = 0
   do np = 1, npc_n
      z(:) = z(:) + tmp(:)
      call MPI_Sendrecv_replace(z, n, MPI_REAL8, np_next, 1111, np_prev, 1111, &
                                mpi_comm_cols, mpi_status, mpierr)
   enddo

end subroutine global_gather

subroutine global_product(z, n)

   ! This routine calculates the global product of z.

   implicit none

   integer n
   real*8 z(n)

   real*8 tmp(n)

   if(npc_n==1 .and. np_rows==1) return ! nothing to do

   ! Do an mpi_allreduce over processor rows

   call mpi_allreduce(z, tmp, n, MPI_REAL8, MPI_PROD, mpi_comm_rows, mpierr)

   ! If only 1 processor column, we are done
   if(npc_n==1) then
      z(:) = tmp(:)
      return
   endif

   ! If all processor columns are involved, we can use mpi_allreduce
   if(npc_n==np_cols) then
      call mpi_allreduce(tmp, z, n, MPI_REAL8, MPI_PROD, mpi_comm_cols, mpierr)
      return
   endif

   ! We send all vectors to the first proc, do the product there
   ! and redistribute the result.

   if(my_pcol == npc_0) then
      z(1:n) = tmp(1:n)
      do np = npc_0+1, npc_0+npc_n-1
         call mpi_recv(tmp,n,MPI_REAL8,np,1111,mpi_comm_cols,mpi_status,mpierr)
         z(1:n) = z(1:n)*tmp(1:n)
      enddo
      do np = npc_0+1, npc_0+npc_n-1
         call mpi_send(z,n,MPI_REAL8,np,1111,mpi_comm_cols,mpierr)
      enddo
   else
      call mpi_send(tmp,n,MPI_REAL8,npc_0,1111,mpi_comm_cols,mpierr)
      call mpi_recv(z  ,n,MPI_REAL8,npc_0,1111,mpi_comm_cols,mpi_status,mpierr)
   endif

end subroutine global_product

subroutine check_monotony(n,d,text)

! This is a test routine for checking if the eigenvalues are monotonically increasing.
! It is for debug purposes only, an error should never be triggered!

   implicit none

   integer n
   real*8 d(n)
   character*(*) text

   integer i

   do i=1,n-1
      if(d(i+1)<d(i)) then
         print '(a,a,i8,2g25.17)','Monotony error on ',text,i,d(i),d(i+1)
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
   enddo

end subroutine check_monotony

end subroutine merge_systems

!-------------------------------------------------------------------------------

subroutine v_add_s(v,n,s)
   implicit none
   integer n
   real*8 v(n),s

   v(:) = v(:) + s
end subroutine v_add_s

!-------------------------------------------------------------------------------

subroutine distribute_global_column(g_col, l_col, noff, nlen, my_prow, np_rows, nblk)

   implicit none

   real*8 g_col(nlen), l_col(*)
   integer noff, nlen, my_prow, np_rows, nblk

   integer nbs, nbe, jb, g_off, l_off, js, je

   nbs = noff/(nblk*np_rows)
   nbe = (noff+nlen-1)/(nblk*np_rows)

   do jb = nbs, nbe

      g_off = jb*nblk*np_rows + nblk*my_prow
      l_off = jb*nblk

      js = MAX(noff+1-g_off,1)
      je = MIN(noff+nlen-g_off,nblk)

      if(je<js) cycle

      l_col(l_off+js:l_off+je) = g_col(g_off+js-noff:g_off+je-noff)

  enddo

end subroutine distribute_global_column

!-------------------------------------------------------------------------------

subroutine solve_secular_equation(n, i, d, z, delta, rho, dlam)

!-------------------------------------------------------------------------------
! This routine solves the secular equation of a symmetric rank 1 modified
! diagonal matrix:
!
!    1. + rho*SUM(z(:)**2/(d(:)-x)) = 0
!
! It does the same as the LAPACK routine DLAED4 but it uses a bisection technique
! which is more robust (it always yields a solution) but also slower
! than the algorithm used in DLAED4.
!
! The same restictions than in DLAED4 hold, namely:
!
!   rho > 0   and   d(i+1) > d(i)
!
! but this routine will not terminate with error if these are not satisfied
! (it will normally converge to a pole in this case).
!
! The output in DELTA(j) is always (D(j) - lambda_I), even for the cases
! N=1 and N=2 which is not compatible with DLAED4.
! Thus this routine shouldn't be used for these cases as a simple replacement
! of DLAED4.
!
! The arguments are the same as in DLAED4 (with the exception of the INFO argument):
!
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension (N)
!         The original eigenvalues.  It is assumed that they are in
!         order, D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (N)
!         DELTA contains (D(j) - lambda_I) in its  j-th component.
!         See remark above about DLAED4 compatibility!
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!-------------------------------------------------------------------------------


   implicit none

   integer n, i
   real*8 d(n), z(n), delta(n), rho, dlam

   integer iter
   real*8 a, b, x, y, dshift

   ! In order to obtain sufficient numerical accuracy we have to shift the problem
   ! either by d(i) or d(i+1), whichever is closer to the solution

   ! Upper and lower bound of the shifted solution interval are a and b


   if(i==n) then

      ! Special case: Last eigenvalue
      ! We shift always by d(n), lower bound is d(n),
      ! upper bound is determined by a guess:

      dshift = d(n)
      delta(:) = d(:) - dshift

      a = 0. ! delta(n)
      b = rho*SUM(z(:)**2) + 1. ! rho*SUM(z(:)**2) is the lower bound for the guess

   else

      ! Other eigenvalues: lower bound is d(i), upper bound is d(i+1)
      ! We check the sign of the function in the midpoint of the interval
      ! in order to determine if eigenvalue is more close to d(i) or d(i+1)

      x = 0.5*(d(i)+d(i+1))
      y = 1. + rho*SUM(z(:)**2/(d(:)-x))

      if(y>0) then
         ! solution is next to d(i)
         dshift = d(i)
      else
         ! solution is next to d(i+1)
         dshift = d(i+1)
      endif

      delta(:) = d(:) - dshift
      a = delta(i)
      b = delta(i+1)

   endif

   ! Bisection:

   do iter=1,200

      ! Interval subdivision

      x = 0.5*(a+b)

      if(x==a .or. x==b) exit   ! No further interval subdivisions possible
      if(abs(x) < 1.d-200) exit ! x next to pole

      ! evaluate value at x

      y = 1. + rho*SUM(z(:)**2/(delta(:)-x))

      if(y==0) then
         ! found exact solution
         exit
      elseif(y>0) then
         b = x
      else
         a = x
      endif

   enddo

   ! Solution:

   dlam = x + dshift
   delta(:) = delta(:) - x

end subroutine solve_secular_equation

!-------------------------------------------------------------------------------

integer function local_index(idx, my_proc, num_procs, nblk, iflag)

!-------------------------------------------------------------------------------
!  local_index: returns the local index for a given global index
!               If the global index has no local index on the
!               processor my_proc behaviour is defined by iflag
!
!  Parameters
!
!  idx         Global index
!
!  my_proc     Processor row/column for which to calculate the local index
!
!  num_procs   Total number of processors along row/column
!
!  nblk        Blocksize
!
!  iflag       Controls the behaviour if idx is not on local processor
!              iflag< 0 : Return last local index before that row/col
!              iflag==0 : Return 0
!              iflag> 0 : Return next local index after that row/col
!-------------------------------------------------------------------------------

   implicit none

   integer idx, my_proc, num_procs, nblk, iflag

   integer iblk

   iblk = (idx-1)/nblk  ! global block number, 0 based

   if(mod(iblk,num_procs) == my_proc) then

      ! block is local, always return local row/col number

      local_index = (iblk/num_procs)*nblk + mod(idx-1,nblk) + 1

   else

      ! non local block

      if(iflag == 0) then

         local_index = 0

      else

         local_index = (iblk/num_procs)*nblk

         if(mod(iblk,num_procs) > my_proc) local_index = local_index + nblk

         if(iflag>0) local_index = local_index + 1
      endif
   endif

end function local_index

!-------------------------------------------------------------------------------

subroutine cholesky_real_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  cholesky_real: Cholesky factorization of a real symmetric matrix
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be factorized.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              On return, the upper triangle contains the Cholesky factor
!              and the lower triangle is set to 0.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
   integer n, nc, i, info
   integer lcs, lce, lrs, lre
   integer tile_size, l_rows_tile, l_cols_tile

   real*8, allocatable:: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   allocate(tmp1(nblk*nblk))
   allocate(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   allocate(tmatr(l_rows,nblk))
   allocate(tmatc(l_cols,nblk))
   tmatr = 0
   tmatc = 0


   do n = 1, na, nblk

      ! Calculate first local row and column of the still remaining matrix
      ! on the local processor

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

      if(n+nblk > na) then

         ! This is the last step, just do a Cholesky-Factorization
         ! of the remaining block

         if(my_prow==prow(n) .and. my_pcol==pcol(n)) then

            call dpotrf('U',na-n+1,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in dpotrf"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

         endif

         exit ! Loop

      endif


      if(my_prow==prow(n)) then

         if(my_pcol==pcol(n)) then

            ! The process owning the upper left remaining block does the
            ! Cholesky-Factorization of this block

            call dpotrf('U',nblk,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in dpotrf"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

            nc = 0
            do i=1,nblk
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            enddo
         endif

         call MPI_Bcast(tmp1,nblk*(nblk+1)/2,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         do i=1,nblk
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         enddo

         if(l_cols-l_colx+1>0) &
            call dtrsm('L','U','T','N',nblk,l_cols-l_colx+1,1.d0,tmp2,ubound(tmp2,1),a(l_row1,l_colx),lda)

      endif

      do i=1,nblk

         if(my_prow==prow(n)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
         if(l_cols-l_colx+1>0) &
            call MPI_Bcast(tmatc(l_colx,i),l_cols-l_colx+1,MPI_REAL8,prow(n),mpi_comm_rows,mpierr)

      enddo

      call elpa_transpose_vectors  (tmatc, ubound(tmatc,1), mpi_comm_cols, &
                                    tmatr, ubound(tmatr,1), mpi_comm_rows, &
                                    n, na, nblk, nblk)

      do i=0,(na-1)/tile_size
         lcs = max(l_colx,i*l_cols_tile+1)
         lce = min(l_cols,(i+1)*l_cols_tile)
         lrs = l_rowx
         lre = min(l_rows,(i+1)*l_rows_tile)
         if(lce<lcs .or. lre<lrs) cycle
         call DGEMM('N','T',lre-lrs+1,lce-lcs+1,nblk,-1.d0, &
                    tmatr(lrs,1),ubound(tmatr,1),tmatc(lcs,1),ubound(tmatc,1), &
                    1.d0,a(lrs,lcs),lda)
      enddo

   enddo

   deallocate(tmp1, tmp2, tmatr, tmatc)

   ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

   do i=1,na
      if(my_pcol==pcol(i)) then
         ! column i is on local processor
         l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
         l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
         a(l_row1:l_rows,l_col1) = 0
      endif
   enddo

end subroutine cholesky_real_2013

!-------------------------------------------------------------------------------

subroutine invert_trm_real_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  invert_trm_real: Inverts a upper triangular matrix
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be inverted.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              The lower triangle is not referenced.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
   integer n, nc, i, info, ns, nb

   real*8, allocatable:: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   allocate(tmp1(nblk*nblk))
   allocate(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   allocate(tmat1(l_rows,nblk))
   allocate(tmat2(nblk,l_cols))
   tmat1 = 0
   tmat2 = 0


   ns = ((na-1)/nblk)*nblk + 1

   do n = ns,1,-nblk

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      nb = nblk
      if(na-n+1 < nblk) nb = na-n+1

      l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)


      if(my_prow==prow(n)) then

         if(my_pcol==pcol(n)) then

            call DTRTRI('U','N',nb,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in DTRTRI"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

            nc = 0
            do i=1,nb
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            enddo
         endif

         call MPI_Bcast(tmp1,nb*(nb+1)/2,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         do i=1,nb
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         enddo

         if(l_cols-l_colx+1>0) &
            call DTRMM('L','U','N','N',nb,l_cols-l_colx+1,1.d0,tmp2,ubound(tmp2,1),a(l_row1,l_colx),lda)

         if(l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
         if(my_pcol==pcol(n)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      endif

      if(l_row1>1) then
         if(my_pcol==pcol(n)) then
            tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
            a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
         endif

         do i=1,nb
            call MPI_Bcast(tmat1(1,i),l_row1-1,MPI_REAL8,pcol(n),mpi_comm_cols,mpierr)
         enddo
      endif

      if(l_cols-l_col1+1>0) &
         call MPI_Bcast(tmat2(1,l_col1),(l_cols-l_col1+1)*nblk,MPI_REAL8,prow(n),mpi_comm_rows,mpierr)

      if(l_row1>1 .and. l_cols-l_col1+1>0) &
         call dgemm('N','N',l_row1-1,l_cols-l_col1+1,nb, -1.d0, &
                    tmat1,ubound(tmat1,1),tmat2(1,l_col1),ubound(tmat2,1), &
                    1.d0, a(1,l_col1),lda)

   enddo

   deallocate(tmp1, tmp2, tmat1, tmat2)

end subroutine invert_trm_real_2013

!-------------------------------------------------------------------------------

subroutine cholesky_complex_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  cholesky_complex: Cholesky factorization of a complex hermitian matrix
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be factorized.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              On return, the upper triangle contains the Cholesky factor
!              and the lower triangle is set to 0.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
   integer n, nc, i, info
   integer lcs, lce, lrs, lre
   integer tile_size, l_rows_tile, l_cols_tile

   complex*16, allocatable:: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
   tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   allocate(tmp1(nblk*nblk))
   allocate(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   allocate(tmatr(l_rows,nblk))
   allocate(tmatc(l_cols,nblk))
   tmatr = 0
   tmatc = 0


   do n = 1, na, nblk

      ! Calculate first local row and column of the still remaining matrix
      ! on the local processor

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

      if(n+nblk > na) then

         ! This is the last step, just do a Cholesky-Factorization
         ! of the remaining block

         if(my_prow==prow(n) .and. my_pcol==pcol(n)) then

            call zpotrf('U',na-n+1,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in zpotrf"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

         endif

         exit ! Loop

      endif


      if(my_prow==prow(n)) then

         if(my_pcol==pcol(n)) then

            ! The process owning the upper left remaining block does the
            ! Cholesky-Factorization of this block

            call zpotrf('U',nblk,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in zpotrf"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

            nc = 0
            do i=1,nblk
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            enddo
         endif

         call MPI_Bcast(tmp1,nblk*(nblk+1)/2,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         do i=1,nblk
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         enddo

         if(l_cols-l_colx+1>0) &
            call ztrsm('L','U','C','N',nblk,l_cols-l_colx+1,(1.d0,0.d0),tmp2,ubound(tmp2,1),a(l_row1,l_colx),lda)

      endif

      do i=1,nblk

         if(my_prow==prow(n)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))
         if(l_cols-l_colx+1>0) &
            call MPI_Bcast(tmatc(l_colx,i),l_cols-l_colx+1,MPI_DOUBLE_COMPLEX,prow(n),mpi_comm_rows,mpierr)

      enddo

      call elpa_transpose_vectors  (tmatc, 2*ubound(tmatc,1), mpi_comm_cols, &
                                    tmatr, 2*ubound(tmatr,1), mpi_comm_rows, &
                                    2*n-1, 2*na, nblk, 2*nblk)


      do i=0,(na-1)/tile_size
         lcs = max(l_colx,i*l_cols_tile+1)
         lce = min(l_cols,(i+1)*l_cols_tile)
         lrs = l_rowx
         lre = min(l_rows,(i+1)*l_rows_tile)
         if(lce<lcs .or. lre<lrs) cycle
         call ZGEMM('N','C',lre-lrs+1,lce-lcs+1,nblk,(-1.d0,0.d0), &
                    tmatr(lrs,1),ubound(tmatr,1),tmatc(lcs,1),ubound(tmatc,1), &
                    (1.d0,0.d0),a(lrs,lcs),lda)
      enddo

   enddo

   deallocate(tmp1, tmp2, tmatr, tmatc)

   ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

   do i=1,na
      if(my_pcol==pcol(i)) then
         ! column i is on local processor
         l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
         l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
         a(l_row1:l_rows,l_col1) = 0
      endif
   enddo

end subroutine cholesky_complex_2013

!-------------------------------------------------------------------------------

subroutine invert_trm_complex_2013(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols)

!-------------------------------------------------------------------------------
!  invert_trm_complex: Inverts a upper triangular matrix
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be inverted.
!              Distribution is like in Scalapack.
!              Only upper triangle is needs to be set.
!              The lower triangle is not referenced.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

   implicit none

   integer na, lda, nblk, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
   integer n, nc, i, info, ns, nb

   complex*16, allocatable:: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)

   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

   allocate(tmp1(nblk*nblk))
   allocate(tmp2(nblk,nblk))
   tmp1 = 0
   tmp2 = 0

   allocate(tmat1(l_rows,nblk))
   allocate(tmat2(nblk,l_cols))
   tmat1 = 0
   tmat2 = 0


   ns = ((na-1)/nblk)*nblk + 1

   do n = ns,1,-nblk

      l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
      l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

      nb = nblk
      if(na-n+1 < nblk) nb = na-n+1

      l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
      l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)


      if(my_prow==prow(n)) then

         if(my_pcol==pcol(n)) then

            call ZTRTRI('U','N',nb,a(l_row1,l_col1),lda,info)
            if(info/=0) then
               print *,"Error in ZTRTRI"
               call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
            endif

            nc = 0
            do i=1,nb
               tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
               nc = nc+i
            enddo
         endif

         call MPI_Bcast(tmp1,nb*(nb+1)/2,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)

         nc = 0
         do i=1,nb
            tmp2(1:i,i) = tmp1(nc+1:nc+i)
            nc = nc+i
         enddo

         if(l_cols-l_colx+1>0) &
            call ZTRMM('L','U','N','N',nb,l_cols-l_colx+1,(1.d0,0.d0),tmp2,ubound(tmp2,1),a(l_row1,l_colx),lda)

         if(l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
         if(my_pcol==pcol(n)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      endif

      if(l_row1>1) then
         if(my_pcol==pcol(n)) then
            tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
            a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
         endif

         do i=1,nb
            call MPI_Bcast(tmat1(1,i),l_row1-1,MPI_DOUBLE_COMPLEX,pcol(n),mpi_comm_cols,mpierr)
         enddo
      endif

      if(l_cols-l_col1+1>0) &
         call MPI_Bcast(tmat2(1,l_col1),(l_cols-l_col1+1)*nblk,MPI_DOUBLE_COMPLEX,prow(n),mpi_comm_rows,mpierr)

      if(l_row1>1 .and. l_cols-l_col1+1>0) &
         call ZGEMM('N','N',l_row1-1,l_cols-l_col1+1,nb, (-1.d0,0.d0), &
                    tmat1,ubound(tmat1,1),tmat2(1,l_col1),ubound(tmat2,1), &
                    (1.d0,0.d0), a(1,l_col1),lda)

   enddo

   deallocate(tmp1, tmp2, tmat1, tmat2)

end subroutine invert_trm_complex_2013

! --------------------------------------------------------------------------------------------------

integer function least_common_multiple(a, b)

   ! Returns the least common multiple of a and b
   ! There may be more efficient ways to do this, we use the most simple approach

   implicit none
   integer, intent(in) :: a, b

   do least_common_multiple = a, a*(b-1), a
      if(mod(least_common_multiple,b)==0) exit
   enddo
   ! if the loop is left regularly, least_common_multiple = a*b

end function

! --------------------------------------------------------------------------------------------------

subroutine hh_transform_real(alpha, xnorm_sq, xf, tau)

   ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
   ! and returns the factor xf by which x has to be scaled.
   ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
   ! since this would be expensive for the parallel implementation.

   implicit none
   real*8, intent(inout) :: alpha
   real*8, intent(in)    :: xnorm_sq
   real*8, intent(out)   :: xf, tau

   real*8 BETA

   if( XNORM_SQ==0. ) then

      if( ALPHA>=0. ) then
         TAU = 0.
      else
         TAU = 2.
         ALPHA = -ALPHA
      endif
      XF = 0.

   else

      BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
      ALPHA = ALPHA + BETA
      IF( BETA<0 ) THEN
         BETA = -BETA
         TAU = -ALPHA / BETA
      ELSE
         ALPHA = XNORM_SQ / ALPHA
         TAU = ALPHA / BETA
         ALPHA = -ALPHA
      END IF
      XF = 1./ALPHA
      ALPHA = BETA

   endif

end subroutine


! --------------------------------------------------------------------------------------------------

subroutine hh_transform_complex(alpha, xnorm_sq, xf, tau)

   ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
   ! and returns the factor xf by which x has to be scaled.
   ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
   ! since this would be expensive for the parallel implementation.

   implicit none
   complex*16, intent(inout) :: alpha
   real*8, intent(in)        :: xnorm_sq
   complex*16, intent(out)   :: xf, tau

   real*8 ALPHR, ALPHI, BETA

   ALPHR = DBLE( ALPHA )
   ALPHI = DIMAG( ALPHA )

   if( XNORM_SQ==0. .AND. ALPHI==0. ) then

      if( ALPHR>=0. ) then
         TAU = 0.
      else
         TAU = 2.
         ALPHA = -ALPHA
      endif
      XF = 0.

   else

      BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
      ALPHA = ALPHA + BETA
      IF( BETA<0 ) THEN
         BETA = -BETA
         TAU = -ALPHA / BETA
      ELSE
         ALPHR = ALPHI * (ALPHI/DBLE( ALPHA ))
         ALPHR = ALPHR + XNORM_SQ/DBLE( ALPHA )
         TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
         ALPHA = DCMPLX( -ALPHR, ALPHI )
      END IF
      XF = 1./ALPHA
      ALPHA = BETA

   endif

end subroutine

! --------------------------------------------------------------------------------------------------

end module ELPA1_2013

! --------------------------------------------------------------------------------------------------
! Please note that the following routines are outside of the module ELPA1
! so that they can be used with real or complex data
! --------------------------------------------------------------------------------------------------

subroutine elpa_transpose_vectors(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvs,nvr,nvc,nblk)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

   use ELPA1_2013 ! for least_common_multiple

   implicit none

   include 'mpif.h'

   integer, intent(in)   :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
   real*8, intent(in)    :: vmat_s(ld_s,nvc)
   real*8, intent(inout) :: vmat_t(ld_t,nvc)

   real*8, allocatable :: aux(:)
   integer myps, mypt, nps, npt
   integer n, lc, k, i, ips, ipt, ns, nl, mpierr
   integer lcm_s_t, nblks_tot, nblks_comm, nblks_skip

   call mpi_comm_rank(comm_s,myps,mpierr)
   call mpi_comm_size(comm_s,nps ,mpierr)
   call mpi_comm_rank(comm_t,mypt,mpierr)
   call mpi_comm_size(comm_t,npt ,mpierr)

   ! The basic idea of this routine is that for every block (in the block cyclic
   ! distribution), the processor within comm_t which owns the diagonal
   ! broadcasts its values of vmat_s to all processors within comm_t.
   ! Of course this has not to be done for every block separately, since
   ! the communictation pattern repeats in the global matrix after
   ! the least common multiple of (nps,npt) blocks

   lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

   nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

   ! Get the number of blocks to be skipped at the begin.
   ! This must be a multiple of lcm_s_t (else it is getting complicated),
   ! thus some elements before nvs will be accessed/set.

   nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

   allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))

   do n = 0, lcm_s_t-1

      ips = mod(n,nps)
      ipt = mod(n,npt)

      if(mypt == ipt) then

         nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
         if(nblks_comm==0) cycle

         if(myps == ips) then
            k = 0
            do lc=1,nvc
               do i = nblks_skip+n, nblks_tot-1, lcm_s_t
                  ns = (i/nps)*nblk ! local start of block i
                  nl = min(nvr-i*nblk,nblk) ! length
                  aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
                  k = k+nblk
               enddo
            enddo
         endif

         call MPI_Bcast(aux,nblks_comm*nblk*nvc,MPI_REAL8,ips,comm_s,mpierr)

         k = 0
         do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
               ns = (i/npt)*nblk ! local start of block i
               nl = min(nvr-i*nblk,nblk) ! length
               vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
               k = k+nblk
            enddo
         enddo

      endif

   enddo

   deallocate(aux)

end subroutine

!-------------------------------------------------------------------------------

subroutine elpa_reduce_add_vectors(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvr,nvc,nblk)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

   use ELPA1_2013 ! for least_common_multiple

   implicit none

   include 'mpif.h'

   integer, intent(in)   :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
   real*8, intent(in)    :: vmat_s(ld_s,nvc)
   real*8, intent(inout) :: vmat_t(ld_t,nvc)

   real*8, allocatable :: aux1(:), aux2(:)
   integer myps, mypt, nps, npt
   integer n, lc, k, i, ips, ipt, ns, nl, mpierr
   integer lcm_s_t, nblks_tot

   call mpi_comm_rank(comm_s,myps,mpierr)
   call mpi_comm_size(comm_s,nps ,mpierr)
   call mpi_comm_rank(comm_t,mypt,mpierr)
   call mpi_comm_size(comm_t,npt ,mpierr)

   ! Look to elpa_transpose_vectors for the basic idea!

   ! The communictation pattern repeats in the global matrix after
   ! the least common multiple of (nps,npt) blocks

   lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

   nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

   allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
   allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
   aux1(:) = 0
   aux2(:) = 0

   do n = 0, lcm_s_t-1

      ips = mod(n,nps)
      ipt = mod(n,npt)

      if(myps == ips) then

         k = 0
         do lc=1,nvc
            do i = n, nblks_tot-1, lcm_s_t
               ns = (i/nps)*nblk ! local start of block i
               nl = min(nvr-i*nblk,nblk) ! length
               aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
               k = k+nblk
            enddo
         enddo

         if(k>0) call mpi_reduce(aux1,aux2,k,MPI_REAL8,MPI_SUM,ipt,comm_t,mpierr)

         if(mypt == ipt) then
            k = 0
            do lc=1,nvc
               do i = n, nblks_tot-1, lcm_s_t
                  ns = (i/npt)*nblk ! local start of block i
                  nl = min(nvr-i*nblk,nblk) ! length
                  vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
                  k = k+nblk
               enddo
            enddo
         endif

      endif

   enddo

   deallocate(aux1)
   deallocate(aux2)

end subroutine

!-------------------------------------------------------------------------------
      SUBROUTINE DLAED4_OWN( N, I, D, Z, DELTA, RHO, DLAM, INFO )
!
!  -- LAPACK routine (version 3.2) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     FHI-aims (VB): We are here providing a verbatim copy of dlaed* from Lapack,
!                    except that we compile it ourselves. In elpa1.f90, we find that
!                    some vendor-compiled Lapack versions can crash in DLAED* with
!                    division-by-zeroes. Our only feasible way to track down these
!                    crashes and avoid them is to include the routines here and
!                    compile them ourselves. Credit for those routines goes to the 
!                    original authors, copyright rests with them, and we believe we
!                    are consistent with their original license (netlib).
!
!     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      DOUBLE PRECISION   DLAM, RHO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the I-th updated eigenvalue of a symmetric
!  rank-one modification to a diagonal matrix whose elements are
!  given in the array d, and that
!
!             D(i) < D(j)  for  i < j
!
!  and that RHO > 0.  This is arranged by the calling routine, and is
!  no loss in generality.  The rank-one modified system is thus
!
!             diag( D )  +  RHO *  Z * Z_transpose.
!
!  where we assume the Euclidean norm of Z is 1.
!
!  The method consists of approximating the rational functions in the
!  secular equation by simpler interpolating rational functions.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension (N)
!         The original eigenvalues.  It is assumed that they are in
!         order, D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (N)
!         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th
!         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5
!         for detail. The vector DELTA contains the information necessary
!         to construct the eigenvectors by DLAED3 and DLAED9.
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!
!  INFO   (output) INTEGER
!         = 0:  successful exit
!         > 0:  if INFO = 1, the updating process failed.
!
!  Internal Parameters
!  ===================
!
!  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!  whether D(i) or D(i+1) is treated as the origin.
!
!            ORGATI = .true.    origin at i
!            ORGATI = .false.   origin at i+1
!
!   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!   if we are working with THREE poles!
!
!   MAXIT is the maximum number of iterations allowed for each
!   eigenvalue.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0, &
                         TEN = 10.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW, &
                         EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI, &
                         RHOINV, TAU, TEMP, TEMP1, W
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZZ( 3 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAED5, DLAED6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
      INFO = 0
      IF( N.EQ.1 ) THEN
!
!         Presumably, I=1 upon entry
!
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL DLAED5_OWN( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      END IF
!
!     Compute machine epsilon
!
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
!
!     The case I = N
!
      IF( I.EQ.N ) THEN
!
!        Initialize some basic variables
!
         II = N - 1
         NITER = 1
!
!        Calculate initial guess
!
         MIDPT = RHO / TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
   10    CONTINUE
!
         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
!
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) + &
             Z( N )*Z( N ) / DELTA( N )
!
         IF( W.LE.ZERO ) THEN
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) + &
                   Z( N )*Z( N ) / RHO
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               IF( A.LT.ZERO ) THEN
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
            END IF
!
!           It can be proved that
!               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
!
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
!
!           It can be proved that
!               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
!
            DLTLB = ZERO
            DLTUB = MIDPT
         END IF
!
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
   30    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                  ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!        Calculate the new step
!
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W - &
             DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO ) &
            C = ABS( C )
         IF( C.EQ.ZERO ) THEN
!          ETA = B/A
!           ETA = RHO - TAU
            ETA = DLTUB - TAU
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF( W*ETA.GT.ZERO ) &
            ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
!
         TAU = TAU + ETA
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                  ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 90 NITER = ITER, MAXIT
!
!           Test for convergence
!
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               DLAM = D( I ) + TAU
               GO TO 250
            END IF
!
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
!
!           Calculate the new step
!
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W - &
                DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF( W*ETA.GT.ZERO ) &
               ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
!
            TAU = TAU + ETA
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                     ABS( TAU )*( DPSI+DPHI )
!
            W = RHOINV + PHI + PSI
   90    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
!
!        End for the case I = N
!
      ELSE
!
!        The case for I < N
!
         NITER = 1
         IP1 = I + 1
!
!        Calculate initial guess
!
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
  100    CONTINUE
!
         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
!
         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) + &
             Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
!
         IF( W.GT.ZERO ) THEN
!
!           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
!
!           We choose d(i) as origin.
!
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            IF( A.GT.ZERO ) THEN
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
!
!           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
!
!           We choose d(i+1) as origin.
!
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = -MIDPT
            DLTUB = ZERO
         END IF
!
         IF( ORGATI ) THEN
!            Original code:
!            DO 130 J = 1, N
!               DELTA( J ) = ( D( J )-D( I ) ) - TAU
!  130       CONTINUE
            ! Terrible things happen if the compiler removes the parenthesis
            ! in the above loop for optimization
            DELTA(1:N) = D(1:N)-D(I)
            DELTA(I) = 0 ! prevent optimizations like loop fusion
            DELTA(1:N) = DELTA(1:N) - TAU
         ELSE
!            DO 140 J = 1, N
!               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
!  140       CONTINUE
            ! See above
            DELTA(1:N) = D(1:N)-D(IP1)
            DELTA(IP1) = 0 ! prevent optimizations like loop fusion
            DELTA(1:N) = DELTA(1:N) - TAU
         END IF
         IF( ORGATI ) THEN
            II = I
         ELSE
            II = I + 1
         END IF
         IIM1 = II - 1
         IIP1 = II + 1
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
!
         W = RHOINV + PHI + PSI
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO ) &
               SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO ) &
               SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N ) &
            SWTCH3 = .FALSE.
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                  THREE*ABS( TEMP ) + ABS( TAU )*DW
!
!        Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!        Calculate the new step
!
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            IF( ORGATI ) THEN
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )* &
                   ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* &
                   ( Z( IP1 ) / DELTA( IP1 ) )**2
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W - &
                DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )* &
                         ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )* &
                         ( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
!
!           Interpolation using THREE most relevant poles
!
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - &
                   ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* &
                         ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - &
                   ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* &
                         ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6_OWN( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, &
                         INFO )
            IF( INFO.NE.0 ) &
               GO TO 250
         END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF( W*ETA.GE.ZERO ) &
            ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
!
         PREW = W
!
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  190    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  200    CONTINUE
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                  THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
!
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN ) &
               SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN ) &
               SWTCH = .TRUE.
         END IF
!
         TAU = TAU + ETA
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 240 NITER = ITER, MAXIT
!
!           Test for convergence
!
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               IF( ORGATI ) THEN
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               END IF
               GO TO 250
            END IF
!
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
!
!           Calculate the new step
!
            IF( .NOT.SWTCH3 ) THEN
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DELTA( IP1 )*DW - &
                         ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* &
                         ( Z( IP1 ) / DELTA( IP1 ) )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               END IF
               A = ( DELTA( I )+DELTA( IP1 ) )*W - &
                   DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DELTA( IP1 )* &
                               DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) + &
                               DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI + &
                            DELTA( IP1 )*DELTA( IP1 )*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
!
!              Interpolation using THREE most relevant poles
!
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - &
                         ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* &
                               ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - &
                         ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* &
                               ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               CALL DLAED6_OWN( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, &
                            INFO )
               IF( INFO.NE.0 ) &
                  GO TO 250
            END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF( W*ETA.GE.ZERO ) &
               ETA = -W / DW
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
!
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
!
            TAU = TAU + ETA
            PREW = W
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  220       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  230       CONTINUE
!
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                     THREE*ABS( TEMP ) + ABS( TAU )*DW
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) &
               SWTCH = .NOT.SWTCH
!
  240    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
!
      END IF
!
  250 CONTINUE
!
      RETURN
!
!     End of DLAED4
!
      END

!-------------------------------------------------------------------------------
      SUBROUTINE DLAED5_OWN( I, D, Z, DELTA, RHO, DLAM )
!
!  -- LAPACK routine (version 3.2) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     FHI-aims (VB): We are here providing a verbatim copy of dlaed* from Lapack,
!                    except that we compile it ourselves. In elpa1.f90, we find that
!                    some vendor-compiled Lapack versions can crash in DLAED* with
!                    division-by-zeroes. Our only feasible way to track down these
!                    crashes and avoid them is to include the routines here and
!                    compile them ourselves. Credit for those routines goes to the 
!                    original authors, copyright rests with them, and we believe we
!                    are consistent with their original license (netlib).
!
!     .. Scalar Arguments ..
      INTEGER            I
      DOUBLE PRECISION   DLAM, RHO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( 2 ), DELTA( 2 ), Z( 2 )
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the I-th eigenvalue of a symmetric rank-one
!  modification of a 2-by-2 diagonal matrix
!
!             diag( D )  +  RHO &  Z * transpose(Z) .
!
!  The diagonal elements in the array D are assumed to satisfy
!
!             D(i) < D(j)  for  i < j .
!
!  We also assume RHO > 0 and that the Euclidean norm of the vector
!  Z is one.
!
!  Arguments
!  =========
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  I = 1 or I = 2.
!
!  D      (input) DOUBLE PRECISION array, dimension (2)
!         The original eigenvalues.  We assume D(1) < D(2).
!
!  Z      (input) DOUBLE PRECISION array, dimension (2)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (2)
!         The vector DELTA contains the information necessary
!         to construct the eigenvectors.
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   B, C, DEL, TAU, TEMP, W
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
      DEL = D( 2 ) - D( 1 )
      IF( I.EQ.1 ) THEN
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         IF( W.GT.ZERO ) THEN
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL
!
!           B > ZERO, always
!
            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            IF( B.GT.ZERO ) THEN
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            ELSE
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            END IF
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         END IF
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
!
!     Now I=2
!
         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         ELSE
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN
!
!     End OF DLAED5
!
      END

!-------------------------------------------------------------------------------
      SUBROUTINE DLAED6_OWN( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
!
!  -- LAPACK routine (version 3.2) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     February 2007
!
!     FHI-aims (VB): We are here providing a verbatim copy of dlaed* from Lapack,
!                    except that we compile it ourselves. In elpa1.f90, we find that
!                    some vendor-compiled Lapack versions can crash in DLAED* with
!                    division-by-zeroes. Our only feasible way to track down these
!                    crashes and avoid them is to include the routines here and
!                    compile them ourselves. Credit for those routines goes to the 
!                    original authors, copyright rests with them, and we believe we
!                    are consistent with their original license (netlib).
!
!     .. Scalar Arguments ..
      LOGICAL            ORGATI
      INTEGER            INFO, KNITER
      DOUBLE PRECISION   FINIT, RHO, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( 3 ), Z( 3 )
!     ..
!
!  Purpose
!  =======
!
!  DLAED6 computes the positive or negative root (closest to the origin)
!  of
!                   z(1)        z(2)        z(3)
!  f(x) =   rho + --------- + ---------- + ---------
!                  d(1)-x      d(2)-x      d(3)-x
!
!  It is assumed that
!
!        if ORGATI = .true. the root is between d(2) and d(3);
!        otherwise it is between d(1) and d(2)
!
!  This routine will be called by DLAED4 when necessary. In most cases,
!  the root sought is the smallest in magnitude, though it might not be
!  in some extremely rare situations.
!
!  Arguments
!  =========
!
!  KNITER       (input) INTEGER
!               Refer to DLAED4 for its significance.
!
!  ORGATI       (input) LOGICAL
!               If ORGATI is true, the needed root is between d(2) and
!               d(3); otherwise it is between d(1) and d(2).  See
!               DLAED4 for further details.
!
!  RHO          (input) DOUBLE PRECISION
!               Refer to the equation f(x) above.
!
!  D            (input) DOUBLE PRECISION array, dimension (3)
!               D satisfies d(1) < d(2) < d(3).
!
!  Z            (input) DOUBLE PRECISION array, dimension (3)
!               Each of the elements in z must be positive.
!
!  FINIT        (input) DOUBLE PRECISION
!               The value of f at 0. It is more accurate than the one
!               evaluated inside this routine (if someone wants to do
!               so).
!
!  TAU          (output) DOUBLE PRECISION
!               The root of the equation f(x).
!
!  INFO         (output) INTEGER
!               = 0: successful exit
!               > 0: if INFO = 1, failure to converge
!
!  Further Details
!  ===============
!
!  30/06/99: Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  10/02/03: This version has a few statements commented out for thread
!  safety (machine parameters are computed on each entry). SJH.
!
!  05/10/06: Modified from a new version of Ren-Cang Li, use
!     Gragg-Thornton-Warner cubic convergent scheme for better stability.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 40 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DSCALE( 3 ), ZSCALE( 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SCALE
      INTEGER            I, ITER, NITER
      DOUBLE PRECISION   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F, &
                         FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1, &
                         SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4,  &
                         LBD, UBD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!

      INFO = 0
!
      IF( ORGATI ) THEN
         LBD = D(2)
         UBD = D(3)
      ELSE
         LBD = D(1)
         UBD = D(2)
      END IF
      IF( FINIT .LT. ZERO )THEN
         LBD = ZERO
      ELSE
         UBD = ZERO 
      END IF
!
      NITER = 1
      TAU = ZERO
      IF( KNITER.EQ.2 ) THEN
         IF( ORGATI ) THEN
            TEMP = ( D( 3 )-D( 2 ) ) / TWO
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
         ELSE
            TEMP = ( D( 1 )-D( 2 ) ) / TWO
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
         END IF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            TAU = B / A
         ELSE IF( A.LE.ZERO ) THEN
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( TAU .LT. LBD .OR. TAU .GT. UBD ) &
            TAU = ( LBD+UBD )/TWO
         IF( D(1).EQ.TAU .OR. D(2).EQ.TAU .OR. D(3).EQ.TAU ) THEN
            TAU = ZERO
         ELSE
            TEMP = FINIT + TAU*Z(1)/( D(1)*( D( 1 )-TAU ) ) + &
                           TAU*Z(2)/( D(2)*( D( 2 )-TAU ) ) + &
                           TAU*Z(3)/( D(3)*( D( 3 )-TAU ) )
            IF( TEMP .LE. ZERO )THEN
               LBD = TAU
            ELSE
               UBD = TAU
            END IF
            IF( ABS( FINIT ).LE.ABS( TEMP ) ) &
               TAU = ZERO
         END IF
      END IF
!
!     get machine parameters for possible scaling to avoid overflow
!
!     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
!     SMINV2, EPS are not SAVEd anymore between one call to the
!     others but recomputed at each call
!
      EPS = DLAMCH( 'Epsilon' )
      BASE = DLAMCH( 'Base' )
      SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) / &
               THREE ) )
      SMINV1 = ONE / SMALL1
      SMALL2 = SMALL1*SMALL1
      SMINV2 = SMINV1*SMINV1
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
      IF( ORGATI ) THEN
         TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
         TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      END IF
      SCALE = .FALSE.
      IF( TEMP.LE.SMALL1 ) THEN
         SCALE = .TRUE.
         IF( TEMP.LE.SMALL2 ) THEN
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
            SCLFAC = SMINV2
            SCLINV = SMALL2
         ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
            SCLFAC = SMINV1
            SCLINV = SMALL1
         END IF
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
         DO 10 I = 1, 3
            DSCALE( I ) = D( I )*SCLFAC
            ZSCALE( I ) = Z( I )*SCLFAC
   10    CONTINUE
         TAU = TAU*SCLFAC
         LBD = LBD*SCLFAC
         UBD = UBD*SCLFAC
      ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
         DO 20 I = 1, 3
            DSCALE( I ) = D( I )
            ZSCALE( I ) = Z( I )
   20    CONTINUE
      END IF
!
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
   30 CONTINUE
      F = FINIT + TAU*FC
!
      IF( ABS( F ).LE.ZERO ) &
         GO TO 60
      IF( F .LE. ZERO )THEN
         LBD = TAU
      ELSE
         UBD = TAU
      END IF
!
!        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
!                            scheme
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
      ITER = NITER + 1
!
      DO 50 NITER = ITER, MAXIT
!
         IF( ORGATI ) THEN
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         END IF
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( F*ETA.GE.ZERO ) THEN
            ETA = -F / DF
         END IF
!
         TAU = TAU + ETA
         IF( TAU .LT. LBD .OR. TAU .GT. UBD ) &
            TAU = ( LBD + UBD )/TWO 
!
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO

!  VB: The loop 40 below has a fix in recent lapack that leads to a correct exit if
!      ( DSCALE( I )-TAU ) == ZERO
!      In my view, this is incorrect; a pole of a function is not its root.
!      I therefore exit incorrectly (info = 1) since elpa will catch incorrect exits
!      anyway by using qr instead of dc.
!
!      However: The code further up fixes the same issue for D(I)-TAU by setting TAU to ZERO .
!      That may be the correct fix, but I have not asserted this yet.
!
!      For now, we exit with info=1 .

         DO 40 I = 1, 3
            IF ( ( DSCALE( I )-TAU ) .ne. ZERO ) THEN
              TEMP = ONE / ( DSCALE( I )-TAU )
              TEMP1 = ZSCALE( I )*TEMP
              TEMP2 = TEMP1*TEMP
              TEMP3 = TEMP2*TEMP
              TEMP4 = TEMP1 / DSCALE( I )
              FC = FC + TEMP4
              ERRETM = ERRETM + ABS( TEMP4 )
              DF = DF + TEMP2
              DDF = DDF + TEMP3
            ELSE
              ! here we conclude INCORRECT convergence.
              INFO = 1
              GO TO 60
            END IF
   40    CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) + &
                  ABS( TAU )*DF
         IF( ABS( F ).LE.EPS*ERRETM ) &
            GO TO 60
         IF( F .LE. ZERO )THEN
            LBD = TAU
         ELSE
            UBD = TAU
         END IF
   50 CONTINUE
      INFO = 1
   60 CONTINUE
!
!     Undo scaling
!
      IF( SCALE ) &
         TAU = TAU*SCLINV
      RETURN
!
!     End of DLAED6
!
      END

!-------------------------------------------------------------------------------
