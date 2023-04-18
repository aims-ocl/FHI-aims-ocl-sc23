Module CC_cl_distribution

  Implicit None
  
  ! Parameters discribe distribution of tensors
  Integer , dimension(:) , allocatable :: CC_mem_ii_D, CC_mem_ii_G
  Integer , dimension(:) , allocatable :: CC_mem_ij_D, CC_mem_ij_G
  Integer , dimension(:) , allocatable :: CC_mem_kl_D
  Integer , dimension(:) , allocatable :: CC_mem_aa_D, CC_mem_aa_G
  Integer , dimension(:) , allocatable :: CC_mem_ab_G, CC_mem_cd_D
  Integer , dimension(:) , allocatable :: CC_mem_bas
  Integer , dimension(:,:) , allocatable :: CC_mem_ab_D

  Integer , dimension(:) , allocatable :: CC_index_ii_D, CC_index_ii_G
  Integer , dimension(:) , allocatable :: CC_index_ij_D, CC_index_ij_G
  Integer , dimension(:) , allocatable :: CC_index_kl_D
  Integer , dimension(:) , allocatable :: CC_index_aa_D, CC_index_aa_G
  Integer , dimension(:) , allocatable :: CC_index_ab_G, CC_index_cd_D
  Integer , dimension(:) , allocatable :: CC_index_bas
  Integer , dimension(:,:) , allocatable :: CC_index_ab_D

  ! Memory distribution for perturbative triples
  Integer , dimension(:) , allocatable :: CC_mem_ijk
  Integer , dimension(:) , allocatable :: CC_index_ijk

  Integer , dimension(:) , allocatable :: CC_mem_abc
  Integer , dimension(:) , allocatable :: CC_index_abc

  Integer (kind=8) :: CC_basic_mem, CC_tmp_mem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

Subroutine CC_cl_vec_distr(distr_start,n_distr,n_tsk,CC_distr_mem,CC_distr_index)

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

End Subroutine CC_cl_vec_distr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module CC_cl_distribution

