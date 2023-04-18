Module CC_3d_distribution

  Implicit None
  
! Parameters discribe distribution of vetors

  ! The four index matrices
  ! Intra-domain distribution of key indices
  Integer , dimension(:,:) , allocatable :: CC_mem_ij_D
  Integer , dimension(:,:) , allocatable :: CC_mem_ab_D
  Integer , dimension(:) , allocatable :: CC_mem_ia_D
  Integer , dimension(:) , allocatable :: CC_mem_cd_D
  Integer , dimension(:) , allocatable :: CC_mem_kl_D
  Integer , dimension(:) , allocatable :: CC_mem_ii_D
  Integer , dimension(:) , allocatable :: CC_mem_aa_D
  Integer , dimension(:) , allocatable :: CC_mem_bas

  Integer , dimension(:,:) , allocatable :: CC_index_ij_D
  Integer , dimension(:,:) , allocatable :: CC_index_ab_D
  Integer , dimension(:) , allocatable :: CC_index_ia_D
  Integer , dimension(:) , allocatable :: CC_index_cd_D
  Integer , dimension(:) , allocatable :: CC_index_kl_D
  Integer , dimension(:) , allocatable :: CC_index_ii_D
  Integer , dimension(:) , allocatable :: CC_index_aa_D
  Integer , dimension(:) , allocatable :: CC_index_bas

  ! Inter-domain distribution of key indices
  Integer , dimension(:) , allocatable :: CC_mem_k1
  Integer , dimension(:) , allocatable :: CC_index_k1

  Integer , dimension(:,:) , allocatable :: CC_mem_k2
  Integer , dimension(:,:) , allocatable :: CC_index_k2

  Integer , dimension(:) , allocatable :: CC_mem_k3
  Integer , dimension(:) , allocatable :: CC_index_k3


  ! The auxiliary vectors below are shared by MPI tasks in same domain,
  ! one replication for each domain.

  ! w matrix, integrals (ab|cd) and b matrix are solved independently by each MPI domain. 
  ! Note that the complete b matrix is never saved, instead of that, b matrix
  ! is divided into many pieces, and each piece it used immediately after
  ! being calculated.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

Subroutine CC_3d_vec_distr(distr_start,n_distr,n_tsk,CC_distr_mem,CC_distr_index)

  Implicit None

  Integer :: distr_start,n_distr,i_tmp,j_tmp
  Integer :: n_tsk, i_task
  Integer , dimension(n_tsk) :: CC_distr_mem, CC_distr_index

  i_tmp = Int(n_distr/n_tsk)
  j_tmp = Mod(n_distr,n_tsk)

  do i_task = 1, n_tsk
    if (i_task.le.j_tmp) then
      CC_distr_mem(i_task) = i_tmp + 1
    else
      CC_distr_mem(i_task) = i_tmp
    end if
  end do

  CC_distr_index(1) = distr_start

  do i_task = 2, n_tsk
    CC_distr_index(i_task) = CC_distr_index(i_task - 1) + CC_distr_mem(i_task - 1)
  end do

End Subroutine CC_3d_vec_distr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module CC_3d_distribution

