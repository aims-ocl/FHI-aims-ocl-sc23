!****h* FHI-aims/restart_elsi
!  NAME
!    restart_elsi
!  SYNOPSIS
module restart_elsi
!  PURPOSE
!    This module provides routines for elsi_restart.
  implicit none
!  AUTHOR
!    Victor Yu, Duke University
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Created, May 2019.
!  SOURCE

  private

  public :: elsi_restart_lapack
  public :: elsi_restart_scalapack
  public :: elsi_restart_check

  interface elsi_restart_lapack
    module procedure elsi_restart_lapack_real
    module procedure elsi_restart_lapack_cmplx
  end interface

  interface elsi_restart_scalapack
    module procedure elsi_restart_scalapack_real
    module procedure elsi_restart_scalapack_cmplx
  end interface

!******
contains
!-------------------------------------------------------------------------------
!****s* restart_elsi/elsi_restart_lapack_real
!  NAME
!    elsi_restart_lapack_real
!  SYNOPSIS
subroutine elsi_restart_lapack_real(i_spin,i_kpt,dm)
!  PURPOSE
!    Reads density matrix and overlap matrix from file.
!    Extrapolates density matrix.
!    P_1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
!    P_0 and S_0 are stored, as they are sparse.
!    U_0 and (U_0 P_0 U_0^T) are dense.
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use dimensions, only: n_basis, n_spin, n_periodic
  use elsi_wrapper, only: rwh_r, aims_elsi_read_mat_dim, aims_elsi_read_mat
  use physics, only: hamiltonian, overlap_matrix
  use runtime_choices, only: elsi_extrap_dm, packed_matrix_format, PM_index
!  SOURCE
  implicit none

  integer, intent(in) :: i_spin
  integer, intent(in) :: i_kpt
  real*8, intent(out) :: dm(n_basis,n_basis)

  character*100 :: mat_file
  real*8 :: n_electrons_file
  integer :: n_basis_file
  integer :: mxld_file
  integer :: mxcol_file
  integer :: ierr
  integer :: i
  integer :: j
  integer :: n
  integer :: nblk
  integer :: nwork

  complex*16 :: dummy1(1,n_spin)
  complex*16 :: dummy2(1)

  real*8, allocatable :: ovlp_work(:,:)
  real*8, allocatable :: work(:,:)
  real*8, allocatable :: ham_tri(:,:)
  real*8, allocatable :: ovlp_tri(:)

  write(mat_file,"(A,I2.2,A,I6.6,A)") "D_spin_",i_spin,"_kpt_",i_kpt,".csc"

  call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
       n_basis_file,mxld_file,mxcol_file)

  call aims_elsi_read_mat(rwh_r,trim(mat_file),dm)

  if(elsi_extrap_dm) then
     write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",i_kpt,".csc"

     call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
          n_basis_file,mxld_file,mxcol_file)

     call aims_allocate(ovlp_work,n_basis,n_basis,"ovlp_work")
     call aims_allocate(work,n_basis,n_basis,"work")

     call aims_elsi_read_mat(rwh_r,trim(mat_file),ovlp_work)

     ! Erase lower triangle of S
     do j = 1,n_basis-1
        do i= j+1,n_basis
           ovlp_work(i,j) = 0.d0
        end do
     end do

     ! ovlp_work = U0
     call dpotrf("U",n_basis,ovlp_work,n_basis,ierr)

     nblk = 128

     ! dm = U_0 P_0 U_0^T
     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call dgemm("N","T",i+nwork-1,nwork,i+nwork-1,1.d0,dm,n_basis,&
             ovlp_work(1,i),n_basis,0.d0,work(1,i),n_basis)
     end do

     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call dgemm("N","N",nwork,n_basis-i+1,i+nwork-1,1.d0,ovlp_work(1,i),&
             n_basis,work(1,i),n_basis,0.d0,dm(i,i),n_basis)
     end do

     ! Construct new overlap
     call aims_allocate(ovlp_tri,n_basis*(n_basis+1)/2,"ovlp_tri")

     if(n_periodic > 0 .or. packed_matrix_format == PM_index) then
        call aims_allocate(ham_tri,n_basis*(n_basis+1)/2,n_spin,"ham_tri")

        call construct_hamiltonian_and_ovl(hamiltonian,overlap_matrix,ham_tri,&
             ovlp_tri,dummy1,dummy2,i_kpt)

        call aims_deallocate(ham_tri,"ham_tri")
     else
        ham_tri = hamiltonian
        ovlp_tri = overlap_matrix
     end if

     ! Unpack overlap
     ovlp_work = 0.d0
     n = 0
     do j = 1,n_basis
        do i = 1,j
           n = n+1
           ovlp_work(i,j) = ovlp_tri(n)
        end do
     end do

     call aims_deallocate(ovlp_tri,"ovlp_tri")

     ! Symmetrize P
     do j = 1,n_basis-1
        do i= j+1,n_basis
           dm(i,j) = dm(j,i)
        end do
     end do

     ! ovlp_work = U1
     call dpotrf("U",n_basis,ovlp_work,n_basis,ierr)

     ! ovlp_work = U1^(-1)
     call dtrtri("U","N",n_basis,ovlp_work,n_basis,ierr)

     ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call dgemm("N","T",i+nwork-1,nwork,i+nwork-1,1.d0,dm,n_basis,&
             ovlp_work(1,i),n_basis,0.d0,work(1,i),n_basis)
     end do

     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call dgemm("N","N",nwork,n_basis-i+1,i+nwork-1,1.d0,ovlp_work(1,i),&
             n_basis,work(1,i),n_basis,0.d0,dm(i,i),n_basis)
     end do

     call aims_deallocate(ovlp_work,"ovlp_work")
     call aims_deallocate(work,"work")
  end if

end subroutine elsi_restart_lapack_real
!******
!-------------------------------------------------------------------------------
!****s* restart_elsi/elsi_restart_lapack_cmplx
!  NAME
!    elsi_restart_lapack_cmplx
!  SYNOPSIS
subroutine elsi_restart_lapack_cmplx(i_spin,i_kpt,dm)
!  PURPOSE
!    Reads density matrix and overlap matrix from file.
!    Extrapolates density matrix.
!    P_1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
!    P_0 and S_0 are stored, as they are sparse.
!    U_0 and (U_0 P_0 U_0^T) are dense.
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use dimensions, only: n_basis, n_spin, n_periodic
  use elsi_wrapper, only: rwh_r, aims_elsi_read_mat_dim, aims_elsi_read_mat
  use physics, only: hamiltonian, overlap_matrix
  use runtime_choices, only: elsi_extrap_dm, packed_matrix_format, PM_index
!  SOURCE
  implicit none

  integer, intent(in) :: i_spin
  integer, intent(in) :: i_kpt
  complex*16, intent(out) :: dm(n_basis,n_basis)

  character*100 :: mat_file
  real*8 :: n_electrons_file
  integer :: n_basis_file
  integer :: mxld_file
  integer :: mxcol_file
  integer :: ierr
  integer :: i
  integer :: j
  integer :: n
  integer :: nblk
  integer :: nwork

  real*8 :: dummy1(1,n_spin)
  real*8 :: dummy2(1)

  complex*16, allocatable :: ovlp_work(:,:)
  complex*16, allocatable :: work(:,:)
  complex*16, allocatable :: ham_tri(:,:)
  complex*16, allocatable :: ovlp_tri(:)

  write(mat_file,"(A,I2.2,A,I6.6,A)") "D_spin_",i_spin,"_kpt_",i_kpt,".csc"

  call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
       n_basis_file,mxld_file,mxcol_file)

  call aims_elsi_read_mat(rwh_r,trim(mat_file),dm)

  if(elsi_extrap_dm) then
     write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",i_kpt,".csc"

     call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
          n_basis_file,mxld_file,mxcol_file)

     call aims_allocate(ovlp_work,n_basis,n_basis,"ovlp_work")
     call aims_allocate(work,n_basis,n_basis,"work")

     call aims_elsi_read_mat(rwh_r,trim(mat_file),ovlp_work)

     ! Erase lower triangle of S
     do j = 1,n_basis-1
        do i= j+1,n_basis
           ovlp_work(i,j) = (0.d0,0.d0)
        end do
     end do

     ! ovlp_work = U0
     call zpotrf("U",n_basis,ovlp_work,n_basis,ierr)

     nblk = 128

     ! dm = U_0 P_0 U_0^T
     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call zgemm("N","C",i+nwork-1,nwork,i+nwork-1,(1.d0,0.d0),dm,n_basis,&
             ovlp_work(1,i),n_basis,(0.d0,0.d0),work(1,i),n_basis)
     end do

     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call zgemm("N","N",nwork,n_basis-i+1,i+nwork-1,(1.d0,0.d0),&
             ovlp_work(1,i),n_basis,work(1,i),n_basis,(0.d0,0.d0),dm(i,i),&
             n_basis)
     end do

     ! Construct new overlap
     call aims_allocate(ovlp_tri,n_basis*(n_basis+1)/2,"ovlp_tri")
     call aims_allocate(ham_tri,n_basis*(n_basis+1)/2,n_spin,"ham_tri")

     call construct_hamiltonian_and_ovl(hamiltonian,overlap_matrix,dummy1,&
          dummy2,ham_tri,ovlp_tri,i_kpt)

     call aims_deallocate(ham_tri,"ham_tri")

     ! Unpack overlap
     ovlp_work = (0.d0,0.d0)
     n = 0
     do j = 1,n_basis
        do i = 1,j
           n = n+1
           ovlp_work(i,j) = ovlp_tri(n)
        end do
     end do

     call aims_deallocate(ovlp_tri,"ovlp_tri")

     ! Symmetrize P
     do j = 1,n_basis-1
        do i= j+1,n_basis
           dm(i,j) = conjg(dm(j,i))
        end do
     end do

     ! ovlp_work = U1
     call zpotrf("U",n_basis,ovlp_work,n_basis,ierr)

     ! ovlp_work = U1^(-1)
     call ztrtri("U","N",n_basis,ovlp_work,n_basis,ierr)

     ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call zgemm("N","C",i+nwork-1,nwork,i+nwork-1,(1.d0,0.d0),dm,n_basis,&
             ovlp_work(1,i),n_basis,(0.d0,0.d0),work(1,i),n_basis)
     end do

     do i = 1,n_basis,nblk
        nwork = nblk

        if(i+nwork-1 > n_basis) then
           nwork = n_basis-i+1
        end if

        call zgemm("N","N",nwork,n_basis-i+1,i+nwork-1,(1.d0,0.d0),&
             ovlp_work(1,i),n_basis,work(1,i),n_basis,(0.d0,0.d0),dm(i,i),&
             n_basis)
     end do

     call aims_deallocate(ovlp_work,"ovlp_work")
     call aims_deallocate(work,"work")
  end if

end subroutine elsi_restart_lapack_cmplx
!******
!-------------------------------------------------------------------------------
!****s* restart_elsi/elsi_restart_scalapack_real
!  NAME
!    elsi_restart_scalapack_real
!  SYNOPSIS
subroutine elsi_restart_scalapack_real(i_spin,dm)
!  PURPOSE
!    Reads density matrix and overlap matrix from file.
!    Extrapolates density matrix.
!    P_1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
!    P_0 and S_0 are stored, as they are sparse.
!    U_0 and (U_0 P_0 U_0^T) are dense.
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use dimensions, only: n_basis
  use elpa1_2013, only: cholesky_real_2013, invert_trm_real_2013, &
      mult_at_b_real_2013
  use elsi_wrapper, only: rwh_r, aims_elsi_read_mat_dim, aims_elsi_read_mat
  use runtime_choices, only: elsi_extrap_dm
  use scalapack_wrapper, only: mxld, mxcol, nprow, npcol, my_scalapack_id, &
      my_k_point, mpi_comm_rows, mpi_comm_cols, sc_desc, ovlp, eigenvec, &
      nb, set_full_matrix_real_L_to_U
!  SOURCE
  implicit none

  integer, intent(in) :: i_spin
  real*8, intent(out) :: dm(mxld,mxcol)

  real*8, allocatable :: old_ovlp(:,:)
  real*8, allocatable :: new_ovlp(:,:)

  character*100 :: mat_file
  real*8 :: n_electrons_file
  integer :: n_basis_file
  integer :: mxld_file
  integer :: mxcol_file

  if(my_scalapack_id < npcol*nprow) then
     write(mat_file,"(A,I2.2,A,I6.6,A)") "D_spin_",i_spin,"_kpt_",my_k_point,&
        ".csc"

     call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
          n_basis_file,mxld_file,mxcol_file)

     call aims_elsi_read_mat(rwh_r,trim(mat_file),dm)

     if(elsi_extrap_dm) then
        write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",my_k_point,".csc"

        call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
             n_basis_file,mxld_file,mxcol_file)

        call aims_allocate(old_ovlp,mxld,mxcol,"old_ovlp")
        call aims_allocate(new_ovlp,mxld,mxcol,"new_ovlp")

        new_ovlp = ovlp

        call aims_elsi_read_mat(rwh_r,trim(mat_file),old_ovlp)

        ! old_ovlp = U0
        call cholesky_real_2013(n_basis,old_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! new_ovlp = U1
        call cholesky_real_2013(n_basis,new_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! new_ovlp = U1^(-1)
        call invert_trm_real_2013(n_basis,new_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! eigenvec = U1^(-T)
        call pdtran(n_basis,n_basis,1.d0,new_ovlp,1,1,sc_desc,0.d0,eigenvec,1,&
             1,sc_desc)

        ! new_ovlp = U_1^(-1) U_0
        call mult_at_b_real_2013("L","U",n_basis,n_basis,eigenvec,mxld,&
             old_ovlp,mxld,nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        ! old_ovlp = U_0^T U_1^(-T)
        call pdtran(n_basis,n_basis,1.d0,new_ovlp,1,1,sc_desc,0.d0,old_ovlp,1,&
             1,sc_desc)

        ! new_ovlp = U_1^(-1) U_0 P_0
        call mult_at_b_real_2013("L","U",n_basis,n_basis,old_ovlp,mxld,dm,mxld,&
             nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        ! dm = P_0 U_0^T U_1^(-T)
        call pdtran(n_basis,n_basis,1.d0,new_ovlp,1,1,sc_desc,0.d0,dm,1,1,&
             sc_desc)

        ! new_ovlp = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
        call mult_at_b_real_2013("L","L",n_basis,n_basis,old_ovlp,mxld,dm,mxld,&
             nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        call set_full_matrix_real_L_to_U(new_ovlp)

        ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
        dm = new_ovlp

        call aims_deallocate(old_ovlp,"old_ovlp")
        call aims_deallocate(new_ovlp,"new_ovlp")
     end if
  end if

end subroutine elsi_restart_scalapack_real
!******
!-------------------------------------------------------------------------------
!****s* restart_elsi/elsi_restart_scalapack_cmplx
!  NAME
!    elsi_restart_scalapack_cmplx
!  SYNOPSIS
subroutine elsi_restart_scalapack_cmplx(i_spin,dm)
!  PURPOSE
!    Reads density matrix and overlap matrix from file.
!    Extrapolates density matrix.
!    P_1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
!    P_0 and S_0 are stored, as they are sparse.
!    U_0 and (U_0 P_0 U_0^T) are dense.
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use dimensions, only: n_basis
  use elpa1_2013, only: cholesky_complex_2013, invert_trm_complex_2013, &
      mult_ah_b_complex_2013
  use elsi_wrapper, only: rwh_r, aims_elsi_read_mat_dim, aims_elsi_read_mat
  use runtime_choices, only: elsi_extrap_dm
  use scalapack_wrapper, only: mxld, mxcol, nprow, npcol, my_scalapack_id, &
      my_k_point, mpi_comm_rows, mpi_comm_cols, sc_desc, ovlp_complex, &
      eigenvec_complex, nb, set_full_matrix_complex_L_to_U
!  SOURCE
  implicit none

  integer, intent(in) :: i_spin
  complex*16, intent(out) :: dm(mxld,mxcol)

  complex*16, allocatable :: old_ovlp(:,:)
  complex*16, allocatable :: new_ovlp(:,:)

  character*100 :: mat_file
  real*8 :: n_electrons_file
  integer :: n_basis_file
  integer :: mxld_file
  integer :: mxcol_file

  if(my_scalapack_id < npcol*nprow) then
     write(mat_file,"(A,I2.2,A,I6.6,A)") "D_spin_",i_spin,"_kpt_",my_k_point,&
        ".csc"

     call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
          n_basis_file,mxld_file,mxcol_file)

     call aims_elsi_read_mat(rwh_r,trim(mat_file),dm)

     if(elsi_extrap_dm) then
        write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",my_k_point,".csc"

        call aims_elsi_read_mat_dim(rwh_r,trim(mat_file),n_electrons_file,&
             n_basis_file,mxld_file,mxcol_file)

        call aims_allocate(old_ovlp,mxld,mxcol,"old_ovlp")
        call aims_allocate(new_ovlp,mxld,mxcol,"new_ovlp")

        new_ovlp = ovlp_complex

        call aims_elsi_read_mat(rwh_r,trim(mat_file),old_ovlp)

        ! old_ovlp = U0
        call cholesky_complex_2013(n_basis,old_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! new_ovlp = U1
        call cholesky_complex_2013(n_basis,new_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! new_ovlp = U1^(-1)
        call invert_trm_complex_2013(n_basis,new_ovlp,mxld,nb,mpi_comm_rows,&
             mpi_comm_cols)

        ! eigenvec = U1^(-T)
        call pztranc(n_basis,n_basis,(1.d0,0.d0),new_ovlp,1,1,sc_desc,&
             (0.d0,0.d0),eigenvec_complex,1,1,sc_desc)

        ! new_ovlp = U_1^(-1) U_0
        call mult_ah_b_complex_2013("L","U",n_basis,n_basis,eigenvec_complex,&
             mxld,old_ovlp,mxld,nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        ! old_ovlp = U_0^T U_1^(-T)
        call pztranc(n_basis,n_basis,(1.d0,0.d0),new_ovlp,1,1,sc_desc,&
             (0.d0,0.d0),old_ovlp,1,1,sc_desc)

        ! new_ovlp = U_1^(-1) U_0 P_0
        call mult_ah_b_complex_2013("L","U",n_basis,n_basis,old_ovlp,mxld,dm,&
             mxld,nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        ! dm = P_0 U_0^T U_1^(-T)
        call pztranc(n_basis,n_basis,(1.d0,0.d0),new_ovlp,1,1,sc_desc,&
             (0.d0,0.d0),dm,1,1,sc_desc)

        ! new_ovlp = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
        call mult_ah_b_complex_2013("L","L",n_basis,n_basis,old_ovlp,mxld,dm,&
             mxld,nb,mpi_comm_rows,mpi_comm_cols,new_ovlp,mxld)

        call set_full_matrix_complex_L_to_U(new_ovlp)

        ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
        dm = new_ovlp

        call aims_deallocate(old_ovlp,"old_ovlp")
        call aims_deallocate(new_ovlp,"new_ovlp")
     end if
  end if

end subroutine elsi_restart_scalapack_cmplx
!******
!-------------------------------------------------------------------------------
!****s* restart_elsi/elsi_restart_check
!  NAME
!    elsi_restart_check
!  SYNOPSIS
  subroutine elsi_restart_check()
!  PURPOSE
!    Stop the code if elsi_restart cannot work.
!  USES
  use dimensions, only: n_spin, n_k_points, n_periodic, n_basis
  use mpi_tasks, only: aims_stop_coll, n_tasks
  use localorb_io, only: localorb_info, use_unit
  use runtime_choices, only: elsi_write_dm, elsi_read_dm, elsi_extrap_dm, &
      use_scalapack
!  SOURCE
  implicit none

  logical :: f_ok
  integer :: bad
  integer :: i_kpt
  integer :: i_spin

  character :: aux1
  character :: aux2
  character :: aux3
  character*100 :: mat_file
  character*200 :: msg

  if(elsi_write_dm .or. elsi_read_dm) then
     if(use_scalapack .and. (n_basis*n_k_points < n_tasks)) then
        write(msg,"(X,A)") "***"
        call localorb_info(msg,use_unit,"(A)")

        if(n_periodic > 0) then
           write(aux1,"(I1)") int(log(real(n_tasks,8))/log(10.d0),4)+1
           write(aux2,"(I1)") int(log(real(n_basis,8))/log(10.d0),4)+1
           write(aux3,"(I1)") int(log(real(n_k_points,8))/log(10.d0),4)+1
           write(msg,"(X,A,I"//aux1//",A,I"//aux2//",A,I"//aux3//",A)") &
              "*** Error - you have ",n_tasks," MPI tasks, ",n_basis,&
              " basis functions, and ",n_k_points," k-points."
           call localorb_info(msg,use_unit,"(A)")
           write(msg,"(X,A)") "***"
           call localorb_info(msg,use_unit,"(A)")
           write(msg,"(X,A)") "*** To use 'elsi_restart', the number of"//&
              " basis functions per k-point cannot be larger than the number"//&
              " of MPI tasks."
           call localorb_info(msg,use_unit,"(A)")
        else
           write(aux1,"(I1)") int(log(real(n_tasks,8))/log(10.d0),4)+1
           write(aux2,"(I1)") int(log(real(n_basis,8))/log(10.d0),4)+1
           write(msg,"(X,A,I"//aux1//",A,I"//aux2//",A)") "*** Error - you"//&
              " have ",n_tasks," MPI tasks and ",n_basis," basis functions."
           call localorb_info(msg,use_unit,"(A)")
           write(msg,"(X,A)") "***"
           call localorb_info(msg,use_unit,"(A)")
           write(msg,"(X,A)") "*** To use 'elsi_restart', the number of"//&
              " basis functions cannot be larger than the number of MPI tasks."
           call localorb_info(msg,use_unit,"(A)")
        end if

        write(msg,"(X,A)") "***"
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(A)") "Too many MPI tasks for 'elsi_restart'."
        call aims_stop_coll(msg)
     end if
  end if

  if(elsi_read_dm) then
     bad = 0

     do i_kpt = 1,n_k_points
        do i_spin = 1,n_spin
           write(mat_file,"(A,I2.2,A,I6.6,A)") "D_spin_",i_spin,"_kpt_",i_kpt,&
              ".csc"

           inquire(file=trim(mat_file),exist=f_ok)

           if(.not. f_ok) then
              if(bad == 0) then
                 write(msg,"(X,A)") "***"
                 call localorb_info(msg,use_unit,"(A)")
              end if

              write(msg,"(X,3A)") "*** Error - the expected file '",&
                 trim(mat_file),"' does not exist."
              call localorb_info(msg,use_unit,"(A)")

              bad = bad+1
           end if
        end do
     end do

     if(elsi_extrap_dm) then
        do i_kpt = 1,n_k_points
           write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",i_kpt,".csc"

           inquire(file=trim(mat_file),exist=f_ok)

           if(.not. f_ok) then
              if(bad == 0) then
                 write(msg,"(X,A)") "***"
                 call localorb_info(msg,use_unit,"(A)")
              end if

              write(msg,"(X,3A)") "*** Error - the expected file '",&
                 trim(mat_file),"' does not exist."
              call localorb_info(msg,use_unit,"(A)")

              bad = bad+1
           end if
        end do
     end if

     if(bad > 0) then
        write(msg,"(X,A)") "***"
        call localorb_info(msg,use_unit,"(A)")
        write(aux1,"(I1)") int(log(real(bad,8))/log(10.d0),4)+1
        write(msg,"(X,A,I"//aux1//",A)") "*** ",bad," file(s) expected by"//&
           " 'elsi_restart read' could not be found."
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(X,A)") "***"
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(X,A)") "*** If you have these files available, copy them"//&
           " to the present working directory."
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(X,A)") "*** Otherwise, remove 'elsi_restart read from'"//&
           " your 'control.in'."
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(X,A)") "***"
        call localorb_info(msg,use_unit,"(A)")
        write(msg,"(A)") "The file(s) needed by 'elsi_restart read' could"//&
           " not be found."
        call aims_stop_coll(msg)
     end if
  end if

  end subroutine elsi_restart_check
!******
end module restart_elsi
