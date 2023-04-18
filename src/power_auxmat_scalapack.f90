!****s* FHI-aims/power_auxmat_scalapack
!  NAME
!   power_auxmat_scalapack
!  SYNOPSIS

      subroutine power_auxmat_scalapack &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a square matrix within the auxiliary basis.
!   This is done by diagonalization.  Eigenmodes corresponding to
!   eigenvalues smaller then prodbas_threshold are set to zero.
!   This is the scalapack version
!
!  USES
       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use scalapack_wrapper
       use localorb_io
       implicit none

!  ARGUMENTS
      real*8, intent(INOUT) :: auxmat(n_basbas,n_loc_prodbas)
      real*8, intent(IN)    :: power
      character*(*), intent(IN) :: name
!
!  INPUTS
!    o auxmat -- real array, e.g. bare coulomb matrix within the auxiliary basis
!    o power -- real number, e.g. -0.5d0
!    o name -- name of matrix (for output only)
!  OUTPUTS
!    o auxmat -- real array, e.g. the inverse square root of the bare coulomb matrix
!                within the auxiliary basis
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      real*8, dimension(:,:), allocatable ::  dist_auxmat
      integer :: info
      character(*), parameter :: func = 'power_auxmat_scalapack'

!  start to work

! Copy auxmat to dist_auxmat changing the matrix distribution from 1D to 2D:

      allocate(dist_auxmat(max_row_2d,max_col_2d), stat=info)
      call check_allocation(info, 'dist_auxmat', func)

      call dist_1d_2d(n_basbas, auxmat, ubound(auxmat,1), dist_auxmat, ubound(dist_auxmat,1))

! diagonalize the matrix

      call power_auxmat_scalapack_2d(dist_auxmat, power, name)

! Redistribute back to 1D

      call dist_2d_1d(n_basbas, dist_auxmat, ubound(dist_auxmat,1), auxmat, ubound(auxmat,1))

!  deallocate arrays

       if(allocated(dist_auxmat)) then
         deallocate(dist_auxmat)
       endif

end subroutine power_auxmat_scalapack
! --------------------------------------------------------------------------------------------------
subroutine dist_1d_2d(nrows,a_1d,ld_a_1d,a_2d,ld_a_2d)

   use dimensions
   use prodbas
   use mpi_tasks

   implicit none

   integer, intent(in) :: nrows, ld_a_1d, ld_a_2d
   real*8, intent(in)  :: a_1d(ld_a_1d,*)
   real*8, intent(out) :: a_2d(ld_a_2d,*)

   integer ic, ip, n, l, i, nown_1d, ic_1d, nown_col_2d, ic_2d, nrows_2d, nrecv, numroc
   integer mpi_status(mpi_status_size), mpierr
   real*8 col(nrows)

   ic_1d = 0
   ic_2d = 0
   nrows_2d = numroc(nrows, nb_aux_2d, myprow_aux_2d, 0, nprow_aux_2d)

   do ic=1,n_basbas

      nown_1d = mod((ic-1)/nb_aux,npcol_aux) ! owner of column ic in 1D distribution
      if(nown_1d == myid) ic_1d = ic_1d + 1
      nown_col_2d = mod((ic-1)/nb_aux_2d,npcol_aux_2d)
      if(nown_col_2d == mypcol_aux_2d) ic_2d = ic_2d + 1

      if(nown_1d /= myid) then

         if(nown_col_2d==mypcol_aux_2d) then
            ! Just receive column ic from PE nown_1d
            call mpi_recv(a_2d(1,ic_2d),nrows_2d,MPI_REAL8,nown_1d,ic,mpi_comm_global,mpi_status,mpierr)
         endif

      else

         ! Send column ic to every PE

         do ip = 0,nprow_aux_2d-1

            n = 0
            do i=ip*nb_aux_2d,nrows-1,nb_aux_2d*nprow_aux_2d
               l = min(nb_aux_2d,nrows-i)
               col(n+1:n+l) = a_1d(i+1:i+l,ic_1d)
               n = n+l
            enddo

            if(n/=numroc(nrows, nb_aux_2d, ip, 0, nprow_aux_2d)) stop 'XXXXXXXXXXXXXXXXXX'

            nrecv = global_id(ip,nown_col_2d)

            if(nrecv==myid) then
               a_2d(1:n,ic_2d) = col(1:n)
            else
               call mpi_send(col,n,MPI_REAL8,nrecv,ic,mpi_comm_global,mpierr)
            endif

         enddo
      endif
   enddo

end subroutine dist_1d_2d

! --------------------------------------------------------------------------------------------------
subroutine dist_2d_1d(nrows,a_2d,ld_a_2d,a_1d,ld_a_1d)

   use dimensions
   use prodbas
   use mpi_tasks

   implicit none

   integer, intent(in) :: nrows, ld_a_1d, ld_a_2d
   real*8, intent(in)  :: a_2d(ld_a_2d,*)
   real*8, intent(out) :: a_1d(ld_a_1d,*)

   integer ic, ip, n, l, i, nown_1d, ic_1d, nown_col_2d, ic_2d, nrows_2d, nrecv, numroc
   integer mpi_status(mpi_status_size), mpierr
   real*8 col(nrows)

   ic_1d = 0
   ic_2d = 0
   nrows_2d = numroc(nrows, nb_aux_2d, myprow_aux_2d, 0, nprow_aux_2d)

   do ic=1,n_basbas

      nown_1d = mod((ic-1)/nb_aux,npcol_aux) ! owner of column ic in 1D distribution
      if(nown_1d == myid) ic_1d = ic_1d + 1
      nown_col_2d = mod((ic-1)/nb_aux_2d,npcol_aux_2d)
      if(nown_col_2d == mypcol_aux_2d) ic_2d = ic_2d + 1

      if(nown_1d /= myid) then

         if(nown_col_2d==mypcol_aux_2d) then
            ! Just send column ic 
            call mpi_send(a_2d(1,ic_2d),nrows_2d,MPI_REAL8,nown_1d,ic,mpi_comm_global,mpierr)
         endif

      else

         ! Recv column ic

         do ip = 0,nprow_aux_2d-1

            n = numroc(nrows, nb_aux_2d, ip, 0, nprow_aux_2d)

            nrecv = global_id(ip,nown_col_2d)

            if(nrecv==myid) then
               col(1:n) = a_2d(1:n,ic_2d)
            else
               call mpi_recv(col,n,MPI_REAL8,nrecv,ic,mpi_comm_global,mpi_status,mpierr)
            endif

            n = 0
            do i=ip*nb_aux_2d,nrows-1,nb_aux_2d*nprow_aux_2d
               l = min(nb_aux_2d,nrows-i)
               a_1d(i+1:i+l,ic_1d) = col(n+1:n+l)
               n = n+l
            enddo

            if(n/=numroc(nrows, nb_aux_2d, ip, 0, nprow_aux_2d)) stop 'YYYYYYYYYYYYYYYYYY'

         enddo
      endif
   enddo

end subroutine dist_2d_1d
! --------------------------------------------------------------------------------------------------

!******
! --------------------------------------------------------------------------------------------------
!****s* FHI-aims/power_auxmat_scalapack
!  NAME
!   power_auxmat_scalapack
!  SYNOPSIS

      subroutine power_auxmat_scalapack_2d &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a square matrix within the auxiliary basis.
!   This is done by diagonalization.  Eigenmodes corresponding to
!   eigenvalues smaller then prodbas_threshold are set to zero.
!   This is the scalapack version
!
!   This version works on a 2-D distributed matrix
!
!  USES
       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use scalapack_wrapper
       use elpa2_2013
       use localorb_io
       implicit none

!  ARGUMENTS
      real*8, intent(INOUT) :: auxmat(max_row_2d,max_col_2d)
      real*8, intent(IN)    :: power
      character*(*), intent(IN) :: name
!
!  INPUTS
!    o auxmat -- real array, e.g. bare coulomb matrix within the auxiliary basis
!    o power -- real number, e.g. -0.5d0
!    o name -- name of matrix (for output only)
!  OUTPUTS
!    o auxmat -- real array, e.g. the inverse square root of the bare coulomb matrix
!                within the auxiliary basis
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      real*8, dimension(:), allocatable ::    eigenvalues
      real*8, dimension(:,:), allocatable ::  transform
      real*8  :: ev_sqrt, threshold

!     working array
      integer   n_nonsingular

      integer   sc_desc_2d(dlen_)
      integer   info
      integer   mpierr
      real*8 :: diff, diff_all
      character*150 :: info_str

!  counter
      integer i_basis_1
      integer lc

!  start to work

      if(myid.eq.0 .and. name /= '') then
       write (use_unit,'(2X,4A)') &
        "Diagonalizing the ", trim(name), " matrix and ", &
        "checking for singularities with ScaLapack ."
       write (use_unit,*)
      endif


      call descinit( sc_desc_2d, n_basbas, n_basbas, nb_aux_2d, nb_aux_2d, &
                     0, 0, &
                     my_blacs_ctxt_aux_2d, MAX(1,max_row_2d), info )

      allocate(transform(max_row_2d,max_col_2d))

! averaging the upper and lower triganles
      call pdtran( n_basbas, n_basbas, 1.0d0, auxmat, &
                   1, 1, sc_desc_2d, 0.0d0, transform, &
                   1, 1, sc_desc_2d )
      if (name /= '') then
         diff = maxval(abs(auxmat - transform))
         call get_max_double(diff_all, diff)
         write(info_str, &
         &     "(2X,'Difference of ',A,' matrix to its transposed:',ES12.4)") &
         & trim(name), diff_all
         call localorb_info(info_str)
      end if

      auxmat(:,:) = 0.5d0* (auxmat(:,:) + transform(:,:))

! diagonalizing the matrix

      if (power < 0.d0) then
         threshold = prodbas_threshold
      else
         ! Still need to cut negative EVs because of the sqrt().
         if (prodbas_threshold.gt.0.d0) then 
           threshold = prodbas_threshold
         else 
           threshold = 0.d0
         end if
      end if

      allocate(eigenvalues(n_basbas))

      auxmat(:,:) = - auxmat(:,:)

      call solve_evp_real_2stage_2013(n_basbas, n_basbas, auxmat, &
           ubound(auxmat,1), eigenvalues, transform, ubound(transform,1), &
           nb_aux_2d, mpi_comm_rows_aux_2d, mpi_comm_cols_aux_2d, &
           mpi_comm_global)

      eigenvalues(1:n_basbas) = -eigenvalues(1:n_basbas)

      do i_basis_1 = 1, n_basbas, 1
         if(eigenvalues(i_basis_1) < threshold) exit
      enddo

      n_nonsingular = i_basis_1-1

      if(myid.eq.0 .and. name /= '') then
         write(use_unit,'(2X,3A,ES14.4,A,ES14.4,A)') &
         & "Eigenvalues of the ", trim(name), " matrix range from", &
         & eigenvalues(1), " to", eigenvalues(n_basbas), "."
         write(use_unit,'(2X,A,I8,A,I8,3A)') &
         & "Using ", n_nonsingular, "   eigenvalues out of rank ", &
         & n_basbas, "   ", trim(name), " matrix (auxiliary basis)."
         if (n_nonsingular < n_basbas) then
            write(use_unit,'(2X,A,ES14.4,A,ES14.4,3A)') &
            & "Still using eigenvalue ", eigenvalues(n_nonsingular), &
            & " while cutting ", eigenvalues(n_nonsingular+1), &
            & " in ", trim(name), " matrix."
         end if
         write(use_unit,*)
      end if
      if (myid.eq.0) then
         if (eigenvalues(n_nonsingular).lt. 1.d-5 .and. power < 0.d0) then
            write(use_unit,'(2X,3A)') &
            & "Be careful! The ", trim(name), " matrix may be ill-conditioned."
            write(use_unit,*)
         end if
      endif

      lc = 0 ! local column counter
      do i_basis_1 = 1, n_nonsingular, 1
         if( MOD((i_basis_1-1)/nb_aux_2d,npcol_aux_2d) .eq. mypcol_aux_2d) then
            ! column i is on local processor
            lc = lc+1
            ev_sqrt = sqrt(eigenvalues(i_basis_1))
            transform(:, lc) = &
               transform(:, lc) * ev_sqrt**power
         endif
      enddo

      call pdgemm('N', 'T', n_basbas, n_basbas, &
                  n_nonsingular, 1.0d0, &
                  transform(1,1), &
                  1, 1, sc_desc_2d, &
                  transform(1,1), &
                  1, 1, sc_desc_2d, 0.d0, &
                  auxmat, &
                  1, 1, sc_desc_2d &
                 )

!  deallocate arrays

       if(allocated(eigenvalues)) then
         deallocate(eigenvalues)
       endif
       if(allocated(transform)) then
         deallocate(transform)
       endif
       end subroutine power_auxmat_scalapack_2d

!******
