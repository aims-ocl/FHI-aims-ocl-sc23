!****s* FHI-aims/scatter_auxmat
!  NAME
!    scatter_auxmat
!  SYNOPSIS

subroutine scatter_auxmat(glob_mat, loc_mat, n_rows)

  !  PURPOSE
  !
  !  USES

  use dimensions
  use mpi_tasks
  use prodbas
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: glob_mat(n_rows, n_basbas)
  real*8, intent(OUT) :: loc_mat(n_rows, n_loc_prodbas)
  integer, intent(IN) :: n_rows

  !  INPUTS
  !    o glob_mat -- On task 0, the auxiliary (e.g. Coulomb) matrix to scatter.
  !    o n_rows -- Number of rows
  !                (n_basbas for coulomb_matr, n_basis_pairs for ovlp_3fn)
  !  OUTPUTS
  !    o loc_mat -- The scattered part for this task.
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: id_send, id_recv, i_task
  integer :: i_basis, i_index
  real*8, allocatable :: mpi_buffer(:,:)
  integer :: mpierr, info, tag
  integer :: my_status(MPI_STATUS_SIZE)

  allocate(mpi_buffer(n_rows, n_max_loc_prodbas), stat=info)
  call check_allocation(info, 'scatter_auxmat:mpi_buffer')

  id_send = 0
  do i_task =  n_tasks, 1, -1
     id_recv = i_task - 1
     tag = id_recv
     if (myid .eq. 0) then
        do i_basis = 1, n_loc_prodbas
           i_index = map_prodbas(i_basis, i_task)
           if(i_index.gt.0) then
              mpi_buffer(:, i_basis) = glob_mat(:, i_index)
           end if
        end do
     end if

     if(id_recv .ne. id_send) then
        if(myid .eq. id_send) then
           call MPI_SEND(mpi_buffer, n_rows*n_max_loc_prodbas, &
           & MPI_DOUBLE_PRECISION, id_recv, tag, mpi_comm_global, mpierr)
        endif
        if(myid .eq. id_recv) then
           call MPI_RECV (mpi_buffer, n_rows*n_max_loc_prodbas, &
           MPI_DOUBLE_PRECISION, id_send, tag, mpi_comm_global,my_status,mpierr)
        endif
     endif
  end do  !  i_task
  loc_mat = mpi_buffer(:, 1:n_loc_prodbas)

  deallocate(mpi_buffer)

end subroutine scatter_auxmat
!******
