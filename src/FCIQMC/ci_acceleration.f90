!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Tonghao Shen
!  SOURCE


  Subroutine ci_acc_optwvfn()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optimize the trial wave-function by solving a eigenvalue   !
! equation H * gamma=lambda * S * gamma, the eigenvector     !
! correspoding to the lowest eigenvalue is used to construct !
! the new trial function for next step                       !
!      See: J. A. Pople, R. Seeger, R. Krishnan,             !
!                       Int. J. Quant. Chem. 11, 149         !
!           R. Seeger, R. Krishnan, J. A. Pople,             !
!                       J. Chem. Phys. 68, 2519              !
!           E. R. Davidson, J. Comp. Phys. 17, 87            !
!           C. D. Sherrill, H. F. Schaefer III               !
!                     Adv. Quant. Chem. 34, 143              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

  implicit none

  Integer :: i_ci_acc,j_ci_acc,k_ci_acc,l_ci_acc,errnum,n_ci_acc
  Double precision , dimension(n_dim_ci_acc_matrix,n_dim_ci_acc_matrix) :: &
             ci_acc_U,ci_acc_UT,ci_acc_mtmp,ci_acc_S_egninv,ci_acc_Hbar,   &
             ci_acc_B_inv

!  Double precision , dimension(:,:) , allocatable :: ci_acc_Hbar

  Double precision , dimension(n_dim_ci_acc_matrix) :: ci_acc_S_egnvl,  &
             ci_acc_gammavec,ci_acc_Hbar_egnvl

  Double precision , dimension(:) , allocatable :: work_ci_acc
  Double precision , dimension(1,1) :: ci_acc_Dummy1

  Integer :: Lwork_ci_acc,info_ci_acc

  Double precision , Parameter :: ci_acc_S_threshold = 1.0E-8
  Double precision :: tmp_coeff_ci

  Logical :: flag_ci_acc_S_singular = .false.

  if ((ci_acc_method.eq.2).or.(ci_acc_method.eq.3)) then

    Call ci_acc_FormHmat()
!    if (myid.eq.0) then
!      write(use_unit,*) 'CI acceleration Form H mat'
!      do i_ci_acc=1,n_dim_ci_acc_matrix
!        write(use_unit,*) (ci_acc_H_matrix(i_ci_acc,j_ci_acc),j_ci_acc=1,n_dim_ci_acc_matrix)
!      end do
!    end if

    n_ci_acc = n_dim_ci_acc_matrix
    ! Solve the eigenvalue problem Hbar * GAMMAprime = lambda * GAMMAprime
    do i_ci_acc = 1, n_ci_acc
      do j_ci_acc = 1, n_ci_acc
        ci_acc_Hbar(i_ci_acc,j_ci_acc) = ci_acc_H_matrix(i_ci_acc,j_ci_acc)
      end do
    end do

    Lwork_ci_acc = n_ci_acc * 3
    Allocate(work_ci_acc(Lwork_ci_acc),stat=errnum)
    Call check_allocation(errnum, 'Work array in CI acceleration')

    Call DSYEV('V','U',n_ci_acc,ci_acc_Hbar,n_ci_acc,ci_acc_Hbar_egnvl, &
               work_ci_acc,Lwork_ci_acc,info_ci_acc)
    Deallocate(work_ci_acc)

!    if (myid.eq.1) then
!      write(use_unit,*) (ci_acc_Hbar_egnvl(j_ci_acc),j_ci_acc=1,n_ci_acc)
!    end if

    ! Construct the folded integral vector according to gamma
    w_0 = 0.0D0
    do j_ci_acc = 1, n_dim_ci_acc_matrix
      w_0 = w_0 + ci_acc_save_w0(j_ci_acc) * ci_acc_Hbar(j_ci_acc,1)
    end do

    w_vect = 0.0D0
    do i_ci_acc = 1, n_configuration_table(myid+1)
      do j_ci_acc = 1, n_dim_ci_acc_matrix
        w_vect(i_ci_acc) = w_vect(i_ci_acc) +  &
                           ci_acc_save_wvect(j_ci_acc,i_ci_acc) * &
                           ci_acc_Hbar(j_ci_acc,1)
      end do
    end do

    ! Construct the optimized wave-function according to gamma
    c_0 = 0.0D0
    do j_ci_acc = 1, n_dim_ci_acc_matrix
      c_0 = c_0 + ci_acc_save_c0(j_ci_acc) * ci_acc_Hbar(j_ci_acc,1)
    end do

    c_vect = 0.0D0
    do i_ci_acc = 1, n_configuration_table(myid+1)
      do j_ci_acc = 1, n_dim_ci_acc_matrix
        c_vect(i_ci_acc) = c_vect(i_ci_acc) +  &
                           ci_acc_save_cvect(j_ci_acc,i_ci_acc) * &
                           ci_acc_Hbar(j_ci_acc,1)
      end do
    end do

!   locate the minimal contribution of basis
    tmp_coeff_ci = 1.0D0
    do i_ci_acc = 1, n_dim_ci_acc_matrix
      if (abs(ci_acc_Hbar(i_ci_acc,1)).lt.tmp_coeff_ci) then
        tmp_coeff_ci = abs(ci_acc_Hbar(i_ci_acc,1))
      end if
    end do

!   th debug
    if (myid.eq.0) then
      write(use_unit,*) 'ci_acc lowest energy',ci_acc_Hbar_egnvl(1)+en_ion_ion
    end if
  end if

  End Subroutine ci_acc_optwvfn


  Subroutine ci_acc_orthogonal(flag_ci_acc_S_singular)

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

  Integer :: n_ci_acc,i_ci_acc,j_ci_acc
  Double precision :: tmp_norm_ci,tmp_dot_ci,ci_acc_norm_fact,tmp_ci_acc_rn
  Double precision , dimension(n_dim_ci_acc_matrix) :: ci_acc_dotpro
  Double precision , Parameter :: ci_acc_S_threshold = 1.0E-4
  Logical :: flag_ci_acc_S_singular

  flag_ci_acc_S_singular = .false.
  ci_acc_dotpro = 0.0D0

  do i_ci_acc = 1, n_dim_ci_acc_matrix
    do j_ci_acc = 1, n_configuration_table(myid+1)
      ci_acc_dotpro(i_ci_acc) = ci_acc_dotpro(i_ci_acc) +           &
                           ci_acc_save_cvect(i_ci_acc,j_ci_acc) *   &
                           c_vect(j_ci_acc)
    end do
    Call sync_real_number(ci_acc_dotpro(i_ci_acc))
    ci_acc_dotpro(i_ci_acc) = ci_acc_dotpro(i_ci_acc) + &
            ci_acc_save_c0(i_ci_acc) * c_0
  end do

  tmp_norm_ci = 1.0D0

  do i_ci_acc = 1, n_dim_ci_acc_matrix
    tmp_norm_ci = tmp_norm_ci - ci_acc_dotpro(i_ci_acc)**2
  end do

!  if (myid.eq.0) then
!    write(use_unit,*) 'ci_acc_dotpro'
!    write(use_unit,*) (ci_acc_dotpro(i_ci_acc),i_ci_acc=1,n_dim_ci_acc_matrix)
!    write(use_unit,*) 'CI acceleration orthogonalization : norm = ',tmp_norm_ci
!  end if

!  n_ci_acc = 0
!  tmp_dot_ci = 0.0D0
!  do i_ci_acc = 2, n_dim_ci_acc_matrix
!    if (abs(ci_acc_dotpro(i_ci_acc)).gt.tmp_dot_ci) then
!      tmp_dot_ci = abs(ci_acc_dotpro(i_ci_acc))
!      n_ci_acc = i_ci_acc
!    end if
!  end do

  if (tmp_norm_ci.lt.ci_acc_S_threshold) then

    flag_ci_acc_S_singular = .true.

    ci_acc_start = ci_acc_start + 1
    if (ci_acc_start.gt.n_dim_ci_acc_matrix) then
      ci_acc_start = 2
    end if

    loc_ci_acc_save = ci_acc_start
    n_ci_acc = ci_acc_start
    tmp_dot_ci = abs(ci_acc_dotpro(n_ci_acc))
    ci_acc_norm_fact = sqrt(1.0D0 / (tmp_norm_ci + tmp_dot_ci**2))

    if (myid.eq.0) then
      write(use_unit,*) 'CI acceleration : S matrix is singular'
    end if

  else

    ci_acc_end = ci_acc_end + 1

    if (ci_acc_end.gt.n_ci_acc_save) then
      ci_acc_start = ci_acc_start + 1
      if (ci_acc_start.gt.n_dim_ci_acc_matrix) then
        ci_acc_start = 2
      end if
      ci_acc_end = n_ci_acc_save
      loc_ci_acc_save = ci_acc_start
      n_ci_acc = ci_acc_start
      tmp_dot_ci = abs(ci_acc_dotpro(n_ci_acc))
      ci_acc_norm_fact = sqrt(1.0D0 / (tmp_norm_ci + tmp_dot_ci**2))
    else
      ci_acc_norm_fact = sqrt(1.0D0 / tmp_norm_ci)
      loc_ci_acc_save = ci_acc_end
      n_ci_acc = 0
    end if

!    if (ci_acc_end.le.n_ci_acc_save) then
!      loc_ci_acc_save = ci_acc_end
!      ci_acc_norm_fact = sqrt(1.0D0 / tmp_norm_ci)
!    else
!      ci_acc_end = n_ci_acc_save
!      tmp_dot_ci = abs(ci_acc_dotpro(loc_ci_acc_save))
!      ci_acc_norm_fact = sqrt(1.0D0 / (tmp_norm_ci + tmp_dot_ci**2))
!    end if

  end if

  if (myid.eq.0) then
    write(use_unit,*) 'loc_ci_acc_save = ',loc_ci_acc_save
  end if

! Project c0
  tmp_ci_acc_rn = c_0
  do j_ci_acc = 1, n_dim_ci_acc_matrix
    if (j_ci_acc.ne.n_ci_acc) then
      tmp_ci_acc_rn = tmp_ci_acc_rn - ci_acc_dotpro(j_ci_acc) &
                                    * ci_acc_save_c0(j_ci_acc)
    end if
  end do
  ci_acc_save_c0(loc_ci_acc_save) = tmp_ci_acc_rn * ci_acc_norm_fact

! Project w0
  tmp_ci_acc_rn = w_0
  do j_ci_acc = 1, n_dim_ci_acc_matrix
    if (j_ci_acc.ne.n_ci_acc) then
      tmp_ci_acc_rn = tmp_ci_acc_rn - ci_acc_dotpro(j_ci_acc) &
                                    * ci_acc_save_w0(j_ci_acc)
    end if
  end do
  ci_acc_save_w0(loc_ci_acc_save) = tmp_ci_acc_rn * ci_acc_norm_fact

! Project coefficints
  do i_ci_acc = 1, n_configuration_table(myid+1)
    tmp_ci_acc_rn = c_vect(i_ci_acc)
    do j_ci_acc = 1, n_dim_ci_acc_matrix
      if (j_ci_acc.ne.n_ci_acc) then
        tmp_ci_acc_rn = tmp_ci_acc_rn - ci_acc_dotpro(j_ci_acc) &
                        * ci_acc_save_cvect(j_ci_acc,i_ci_acc)
      end if
    end do
    ci_acc_save_cvect(loc_ci_acc_save,i_ci_acc) = &
                                 tmp_ci_acc_rn * ci_acc_norm_fact
  end do

! Project integrals
  do i_ci_acc = 1, n_configuration_table(myid+1)
    tmp_ci_acc_rn = w_vect(i_ci_acc)
    do j_ci_acc = 1, n_dim_ci_acc_matrix
      if (j_ci_acc.ne.n_ci_acc) then
        tmp_ci_acc_rn = tmp_ci_acc_rn - ci_acc_dotpro(j_ci_acc) &
                        * ci_acc_save_wvect(j_ci_acc,i_ci_acc)
      end if
    end do
    ci_acc_save_wvect(loc_ci_acc_save,i_ci_acc) = &
                                 tmp_ci_acc_rn * ci_acc_norm_fact
  end do


  End Subroutine ci_acc_orthogonal


  Subroutine ci_acc_FormHmat()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct the Hamiltonian matrix in the subspace  !
! consist of trial wave functions                   !
! H_ij = <psi_1|H|psi_j>                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

  Integer :: i_ci_acc,j_ci_acc,k_ci_acc
  Double precision :: h_ci_acc_tmp

  ! Since H matrix is symmetric, calculating upper triangular part is sufficient
  do i_ci_acc = 1, n_dim_ci_acc_matrix
    do j_ci_acc = i_ci_acc, n_dim_ci_acc_matrix

      h_ci_acc_tmp = 0.0D0

      do k_ci_acc = 1, n_configuration_table(myid+1)
        h_ci_acc_tmp = h_ci_acc_tmp + ci_acc_save_cvect(i_ci_acc,k_ci_acc) &
                      * ((ci_acc_Es(k_ci_acc) + V_00)                      &
                      * ci_acc_save_cvect(j_ci_acc,k_ci_acc)               &
                      + ci_acc_save_wvect(j_ci_acc,k_ci_acc))
      end do

      Call sync_real_number(h_ci_acc_tmp)

      h_ci_acc_tmp = h_ci_acc_tmp + ci_acc_save_c0(i_ci_acc) &
                     * (E_HF * ci_acc_save_c0(j_ci_acc) +    &
                        ci_acc_save_w0(j_ci_acc))

      ! determine H_ij and H_ji
      ci_acc_H_matrix(i_ci_acc,j_ci_acc) = h_ci_acc_tmp

      if (i_ci_acc.ne.j_ci_acc) then
        ci_acc_H_matrix(j_ci_acc,i_ci_acc) = h_ci_acc_tmp
      end if

    end do
  end do

  End Subroutine ci_acc_FormHmat




  Subroutine ci_acc_FormSmat()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct overlap matrix in CI SCF procedure      !
! S(i,j)=<psi_i|psi_j>                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

  Integer :: i_ci_acc,j_ci_acc,k_ci_acc,l_ci_acc
  Double precision :: s_ci_acc_tmp

  do i_ci_acc = 1, n_dim_ci_acc_matrix

!    ci_acc_ovlp_matrix(i_ci_acc,i_ci_acc) = 1.0D0

!    do j_ci_acc = i_ci_acc+1, n_dim_ci_acc_matrix
    do j_ci_acc = i_ci_acc, n_dim_ci_acc_matrix

      s_ci_acc_tmp = 0.0D0

      do k_ci_acc = 1, n_configuration_table(myid+1)
        s_ci_acc_tmp = s_ci_acc_tmp + ci_acc_save_cvect(i_ci_acc,k_ci_acc) &
                                    * ci_acc_save_cvect(j_ci_acc,k_ci_acc)
      end do

      Call sync_real_number(s_ci_acc_tmp)

      s_ci_acc_tmp = s_ci_acc_tmp + ci_acc_save_c0(i_ci_acc) &
                                  * ci_acc_save_c0(j_ci_acc)

      ! determine S_ij and S_ji
      ci_acc_ovlp_matrix(i_ci_acc,j_ci_acc) = s_ci_acc_tmp

      if (i_ci_acc.ne.j_ci_acc) then
        ci_acc_ovlp_matrix(j_ci_acc,i_ci_acc) = s_ci_acc_tmp
      end if

    end do
  end do

  End Subroutine ci_acc_FormSmat


