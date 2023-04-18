Module CC_3d_mem

  Implicit None

  Double complex , dimension(:,:,:) , allocatable :: CC_t_s
! Single excitation coefficients
! CC_t_s(i,a,ki=ka)

  Double complex , dimension(:,:,:) , allocatable :: CC_w_s
! Working vector for single excitations
! CC_w_s(i,a,ki=ka)

  Double complex , dimension(:,:,:) , allocatable :: CC_h_ik, CC_h_ca, CC_h_ck
! Auxiliary vectors

  Double complex , dimension(:,:,:) , allocatable :: CC_g_ca, CC_g_ik
! Auxiliary vectors

  Double complex , dimension(:,:) , allocatable :: CC_DIIS_Bmat
  Double complex , dimension(:) , allocatable :: CC_DIIS_tau
! DIIS vectors

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_t_d
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_t_d_A
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_t_d_nd
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_t_d_d

! Double excitation coefficients

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_w_d_nd
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_w_d_d
! Working vector for double excitations

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_a_aux
! auxiliary vector a

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_j_aux
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_k_aux
! auxiliary vector j and k

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_ten_con
! tensor contraction result

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_aikc
! integral tensor (ai|kc) (for j aux, k_aux and t(i,a))

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_kiac
! integral tensor (ki|ac)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_ackd
! integral tensor (ac|kd) ---- W (for j aux and k aux)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_ackd_nd
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_intl_ackd_d
! integral tensor (ac|kd) for b aux

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_kilc
! integral tensor (ki|lc)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_likc
! integral tensor (li|kc)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_kilc_T
! integral tensor (ki|lc)

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_intl_kilj
! integral tensor (ki|lj) for a aux

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_intl_acbd_nd, CC_intl_acbd_d
! integral tensor (ac|bd) for b aux (optional)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_aibc
! integral tensor (ai|bc)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_aikj
! integral tensor (ai|kj)

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_kcld
! integral tensor (kc|ld) 

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_intl_kcld_A
! integral tensor (kc|ld) ---- A

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: CC_intl_iajb_nd
! integral tensor (ia|jb) ---- nd

  double complex , dimension(:,:,:,:,:,:) , allocatable :: cc_intl_iajb_d
! integral tensor (ia|jb) ---- d

  Double complex , dimension(:,:,:,:,:) , allocatable :: CC_RI_L, CC_RI_R
! RI coefficients

  Double complex , dimension(:,:,:,:) , allocatable :: CC_DIIS_t_s, CC_DIIS_r_s
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: CC_DIIS_t_d, CC_DIIS_t_nd, &
                                                             CC_DIIS_r_d, CC_DIIS_r_nd
! vectors used in DIIS procedure


  Double complex , dimension(:,:,:,:) , allocatable :: CC_t_d_current
! Current coefficients for double excitation (diagonal)

End Module CC_3d_mem

