  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_coulomb_matr_recip
  !  NAME
  !    get_coulomb_matr_recip
  !  SYNOPSIS

  subroutine get_coulomb_matr_recip &
            (coulomb_matr_recip, output_info)

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
    use mpi_tasks, only: check_allocation, myid, n_tasks
    implicit none

    !  ARGUMENTS

       complex*16, intent(OUT) :: coulomb_matr_recip(n_basbas,n_basbas,n_q_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_matr_recip: Coulomb matrix in reciprocal sapce
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8, allocatable :: k_lattvec(:,:) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_matr_recip'

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating the Coulomb matrix in reciprocal space ..."
      call localorb_info(info_str)
    endif

   allocate(k_lattvec(3,n_q_points_task),stat=info)
   call check_allocation(info, 'k_lattvect', func)


    k_lattvec(:,:) = 0.d0
    do i_k_point = 1, n_q_points, 1

        if(myid.eq.mod(i_k_point,n_tasks)) then
         
          i_k_point_local = (i_k_point-1)/n_tasks + 1
          k_lattvec(1:3,i_k_point_local) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))

        endif

    enddo

    call get_qspace_auxmat(n_q_points_task,k_lattvec,coulomb_matr_recip)

    if(allocated(k_lattvec)) then
      deallocate(k_lattvec)
    endif

  end subroutine get_coulomb_matr_recip


  subroutine get_coulomb_matr_blacs(coulomb_matr_blacs, output_info)

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
    implicit none

    !  ARGUMENTS

       complex*16, intent(OUT) :: coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col,n_irkq_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_matr_recip: Coulomb matrix in reciprocal sapce
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8 :: k_lattvec(3,n_irkq_points_task) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_irk_point_local, i_irk_point
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_matr_recip'

    call perfon('gcou')

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating the Coulomb matrix in reciprocal space ..."
      call localorb_info(info_str)
    endif

    k_lattvec(:,:) = 0.d0
    do i_irk_point_local = 1, n_irkq_points_task
       i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
       i_k_point = inv_irk_point_mapping(i_irk_point)

       k_lattvec(1:3,i_irk_point_local) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))
    enddo

    call get_qspace_auxmat(n_irkq_points_task,k_lattvec,coulomb_matr_blacs,lbb_row-1,n_bb_row,lbb_col-1,n_bb_col)
    call perfoff

  end subroutine get_coulomb_matr_blacs

  subroutine get_full_coulomb_matr_blacs(coulomb_matr_blacs, output_info)

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
    implicit none

    !  ARGUMENTS

       complex*16, intent(OUT) :: coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col,n_ks_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_matr_recip: Coulomb matrix in reciprocal sapce
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8 :: k_lattvec(3,n_ks_points_task) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_matr_recip'

    call perfon('gcou')

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating the Coulomb matrix in reciprocal space ..."
      call localorb_info(info_str)
    endif

    k_lattvec(:,:) = 0.d0
    do i_k_point_local = 1, n_ks_points_task
       i_k_point=n_tasks_irkq*(i_k_point_local-1) + myid_irkq + 1

       k_lattvec(1:3,i_k_point_local) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))
    enddo

    call get_qspace_auxmat(n_ks_points_task,k_lattvec,coulomb_matr_blacs,lbb_row-1,n_bb_row,lbb_col-1,n_bb_col)
    call perfoff

  end subroutine get_full_coulomb_matr_blacs

  subroutine get_coulomb_matr_blacs_cpt2(coulomb_matr_blacs, output_info)

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
    use cpt2_blacs
    implicit none

    !  ARGUMENTS

       complex*16, intent(OUT) :: coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col,n_kq_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_matr_recip: Coulomb matrix in reciprocal sapce
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8 :: k_lattvec(3,n_kq_points_task) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_matr_recip'

    call perfon('gcou')

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating the Coulomb matrix in reciprocal space ..."
      call localorb_info(info_str)
    endif

    k_lattvec(:,:) = 0.d0
    do i_k_point_local = 1, n_kq_points_task
       i_k_point=n_tasks_kq_cpt2*(i_k_point_local-1) + myid_kq_cpt2 + 1
       !i_k_point = inv_irk_point_mapping(i_irk_point)

       k_lattvec(1:3,i_k_point_local) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))
    enddo

    call get_qspace_auxmat(n_kq_points_task,k_lattvec,coulomb_matr_blacs,lbb_row-1,n_bb_row,lbb_col-1,n_bb_col)
    call perfoff

  end subroutine get_coulomb_matr_blacs_cpt2

