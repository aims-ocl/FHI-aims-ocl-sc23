  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_coulomb_coeff_blacs
  !  NAME
  !    get_coulomb_coeff_blacs
  !  SYNOPSIS

  subroutine get_coulomb_coeff_blacs &
            (coulomb_coeff_blacs, output_info)

    !  PURPOSE
    !
    !    Fourier transfrom the LVL triple coefficients from real space (on a Bravais lattice)
    !    to reciprocal space.
    !   
    !
    !  USES

    use dimensions
    use prodbas
    use pbc_lists
    use geometry
    use tight_binding_auxmat
    use localorb_io
     use crpa_blacs
    use mpi_tasks
    implicit none

    !  ARGUMENTS

       complex*16, intent(OUT) :: coulomb_coeff_blacs(lbb_row:ubb_row,lbb_col:ubb_col,n_irkq_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_coeff_blacs: the coeffients of 1/q^2 contributions of the Coulomb matrix at q->0
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8  qsq
    real*8, allocatable :: k_lattvec(:,:) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_irk_point, i_irk_point_local
    integer :: i_basis_1, i_basis_2
    integer :: i_basis_fn_1, i_basis_fn_2
    integer :: i_task_gamma
    integer :: info
    integer :: mpierr
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_coeff_blacs'

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating coefficents of the Coulomb matrix at small q vector ..."
      call localorb_info(info_str)
    endif

   allocate(k_lattvec(3,n_irkq_points_task),stat=info)
   call check_allocation(info, 'k_lattvec', func)


    k_lattvec(:,:) = 0.d0
    do i_irk_point_local = 1, n_irkq_points_task
       i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
       i_k_point = inv_irk_point_mapping(i_irk_point)

       k_lattvec(1:3,i_irk_point_local) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))

       if(all(abs(k_point_list(i_k_point,:)).lt.1.e-10)) then
            i_task_gamma = mod(i_k_point, n_tasks)
       endif
    enddo



    call get_qspace_coulomb_coeff(n_irkq_points_task,k_lattvec,coulomb_coeff_blacs,lbb_row-1,n_bb_row,lbb_col-1,n_bb_col)
!    if(n_tasks .gt. 1) then
!      call mpi_bcast(coulomb_coeff_blacs,n_basbas*n_basbas,MPI_COMPLEX16,i_task_gamma,mpi_comm_global, mpierr)
!    endif
    do i_irk_point_local = 1, n_irkq_points_task, 1
        do i_basis_2 = lbb_col, ubb_col, 1
         do i_basis_1 = lbb_row, ubb_row, 1
              if(basbas_l(i_basis_1).ne. 0 .or. basbas_l(i_basis_2).ne. 0) then
                 coulomb_coeff_blacs(i_basis_1,i_basis_2,i_irk_point_local)=(0.d0,0.d0)
              endif
           enddo
         enddo

        qsq=k_lattvec(1,i_irk_point_local)**2 + k_lattvec(2,i_irk_point_local)**2 + k_lattvec(3,i_irk_point_local)**2
        i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
        i_k_point = inv_irk_point_mapping(i_irk_point)

        if(abs(k_point_list(i_k_point,1)).gt.1.e-10 .or. abs(k_point_list(i_k_point,2)).gt.1.e-10 ) cycle
        do i_basis_2=lbb_col, ubb_col, 1
          i_basis_fn_2=basbas_fn(i_basis_2)
          do i_basis_1=lbb_row, ubb_row, 1
             i_basis_fn_1=basbas_fn(i_basis_1)
             if((basbas_l(i_basis_1)+ basbas_l(i_basis_2).le.1) .and. multipole_basbas_fn(i_basis_fn_1) .gt. 1.e-10 &
                         .and. multipole_basbas_fn(i_basis_fn_2).gt.1.e-10 ) then
              write(use_unit,'(I4, f13.8, 4I5, 2f18.8)') i_irk_point_local, sqrt(qsq), i_basis_1, i_basis_2, basbas_l(i_basis_1), basbas_l(i_basis_2), coulomb_coeff_blacs(i_basis_1,i_basis_2,i_irk_point_local)
            endif
          enddo
        enddo
    enddo


    if(allocated(k_lattvec)) then
      deallocate(k_lattvec)
    endif

  end subroutine get_coulomb_coeff_blacs
