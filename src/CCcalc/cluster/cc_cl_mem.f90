Module CC_cl_mem

  Implicit None

  Double precision , dimension(:,:) , allocatable :: CC_t_s
! Single excitation coefficients
! CC_t_s(i,a)

  Double precision , dimension(:,:) , allocatable :: CC_w_s
! Working vector for single excitations

  Double precision , dimension(:,:) , allocatable :: CC_h_ik, CC_h_ca, CC_h_ck
! Auxiliary vectors

  Double precision , dimension(:,:) , allocatable :: CC_g_ca, CC_g_ik
! Auxiliary vectors

  Double precision , dimension(:,:) , allocatable :: CC_RLE_Rmat
  Double precision , dimension(:) , allocatable :: CC_RLE_tau
! RLE matrix and vector  

  Double precision , dimension(:,:) , allocatable :: CC_DIIS_Bmat
  Double precision , dimension(:) , allocatable :: CC_DIIS_tau
! DIIS matrix and vector  

  Double precision , dimension(:,:,:,:) , allocatable :: CC_t_d
  Double precision , dimension(:,:,:,:) , allocatable :: CC_t_d_T,CC_t_d_I
  Double precision , dimension(:,:,:,:) , allocatable :: CC_tau_m, CC_tau_p
! Double excitation coefficients

  Double precision , dimension(:,:,:,:) , allocatable :: CC_w_d
! Working vector for double excitations

  Double precision , dimension(:,:,:,:) , allocatable :: CC_ten_con
! Tensor contraction result

  Double precision , dimension(:,:,:,:) , allocatable :: CC_a_aux
! auxiliary vector a

  Double precision , dimension(:,:,:,:) , allocatable :: CC_j_aux
  Double precision , dimension(:,:,:,:) , allocatable :: CC_k_aux
! auxiliary vector j and k

  Double precision , dimension(:,:,:,:) , allocatable :: CC_intl_iajb, CC_intl_iajb_A

  Double precision , dimension(:,:,:,:) , allocatable :: CC_intl_kiac, CC_intl_likc, &
                                                         CC_intl_kilc, CC_intl_kilc_I, &
                                                         CC_intl_kilj, CC_intl_acbd

  Double precision , dimension(:,:,:,:) , allocatable :: CC_intl_ackd

  Double precision , dimension(:,:,:,:) , allocatable :: CC_intl_dbia,CC_intl_jlkc, &
                                                         CC_intl_iajb_I
! orbital integrals

  Double precision , dimension(:,:,:,:) , allocatable :: CC_RI
  Double precision , dimension(:,:,:) , allocatable :: CC_RI_B
! RI coefficients

  Double precision , dimension(:,:,:,:) , allocatable :: CC_DIIS_t_s, CC_DIIS_r_s
  Double precision , dimension(:,:,:,:) , allocatable :: CC_DIIS_t_d, CC_DIIS_t_nd, &
                                                         CC_DIIS_r_d, CC_DIIS_r_nd
! vectors used in DIIS procedure

  Double precision , dimension(:,:,:) , allocatable :: CC_PT_W_abc

  Double precision , dimension(:,:,:,:) , allocatable :: CC_PT_U

  Double precision , dimension(:,:,:,:) , allocatable :: CC_PT_T_abc
  Double precision , dimension(:,:,:,:,:) , allocatable :: CC_PT_T_send, CC_PT_T_recv

  Double precision , dimension(:,:,:) , allocatable :: CC_PT_jlkc_abc
  Double precision , dimension(:,:,:,:) , allocatable :: CC_PT_jlkc_send, CC_PT_jlkc_recv

  Double precision , dimension(:,:,:) , allocatable :: CC_PT_iajb_abc
  Double precision , dimension(:,:,:,:) , allocatable :: CC_PT_iajb_send, CC_PT_iajb_recv

! vectors used in perturbative T correction

  Double precision , dimension(:,:,:,:) , allocatable :: CC_t_check,CC_intl_bdai_check, &
                                                         CC_intl_ckjl_check, CC_intl_bjck_check
! used for debug

End Module CC_cl_mem

