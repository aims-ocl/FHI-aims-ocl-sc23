subroutine coulomb_4ks(val, cond, ovlp_3ks_bse, coulomb_4ks_vcvc, coulomb_4ks_vccv)
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  implicit none
! arguments
  integer :: val, cond
  real*8, dimension(n_basbas,n_states,n_states,n_spin) :: ovlp_3ks_bse
  real*8, dimension(val, cond, val, cond) :: coulomb_4ks_vcvc
  real*8, dimension(val, cond, cond, val) :: coulomb_4ks_vccv
! local variables
  real*8 :: temp_4ks
  real*8, dimension(:, :), allocatable :: ovlp_3ks_vc, ovlp_3ks_cv, ovlp_4ks_temp
  integer :: v1, v2, c1, c2, counter, col, row, prod_vc

  prod_vc = val * cond
  allocate(ovlp_3ks_vc(n_basbas, prod_vc))
  do row = 1, n_basbas
    do v1 = 1, val
      do c1 = val + 1, n_states
        ovlp_3ks_vc(row, (v1 - 1)*cond + c1-val)  = ovlp_3ks_bse(row, v1, c1, n_spin)
      end do
    end do
  end do
  allocate(ovlp_4ks_temp(prod_vc, prod_vc))
  ovlp_4ks_temp = 0
  call dgemm('T', 'N', prod_vc, prod_vc, n_basbas, 1.0d0,&
             ovlp_3ks_vc, n_basbas, ovlp_3ks_vc, n_basbas, 0.0d0, ovlp_4ks_temp, val * cond)
!  open(unit = 99, file = "coulomb_4ks_vcvc_dgemm")
  coulomb_4ks_vcvc = 0
  do v1 = 1, val
    do c1 = 1, cond
      do v2 = 1, val
        do c2 = 1, cond
          coulomb_4ks_vcvc(v1, c1, v2, c2) = ovlp_4ks_temp((v1 - 1)*cond + c1, (v2 - 1)*cond + c2) 
!          if(abs(coulomb_4ks_vcvc(v1, c1, v2, c2)) > 1e-9) write(99, *) v1, c1, v2, c2, coulomb_4ks_vcvc(v1, c1, v2, c2)
        end do
      end do
    end do
  end do
!  close(99)

  allocate(ovlp_3ks_cv(n_basbas, prod_vc))
  do row = 1, n_basbas
    do c1 = val + 1, n_states
      do v1 = 1, val
        ovlp_3ks_cv(row, (c1 - val - 1) * val + v1)  = ovlp_3ks_bse(row, c1, v1, n_spin)
      end do
    end do
  end do
  ovlp_4ks_temp = 0
  call dgemm('T', 'N', prod_vc, prod_vc, n_basbas, 1.0d0,&
             ovlp_3ks_vc, n_basbas, ovlp_3ks_cv, n_basbas, 0.0d0, ovlp_4ks_temp, val * cond)
!  open(unit = 99, file = "coulomb_4ks_vccv_dgemm")
  coulomb_4ks_vccv = 0
  do v1 = 1, val
    do c1 = 1, cond
      do c2 = 1, cond
        do v2 = 1, val
          coulomb_4ks_vccv(v1, c1, c2, v2) = ovlp_4ks_temp((v1 - 1) * cond + c1, (c2 - 1) * val + v2)
!          if(abs(coulomb_4ks_vccv(v1, c1, c2, v2)) > 1e-9) write(99, *) v1, c1, c2, v2, coulomb_4ks_vccv(v1, c1, c2, v2)
        end do
      end do
    end do
  end do
!  close(99)
!  open(unit = 99, file = "coulomb_4ks_vcvc_dgemm")
!  do v1 = 1, val
!    do c1 = val + 1, n_states
!      do v2 = 1, val
!        do c2 = val + 1, n_states
!          temp_4ks = 0
!          do counter = 1, n_basbas
!            temp_4ks = temp_4ks + ovlp_3ks_bse(counter, v1, c1, n_spin) * ovlp_3ks_bse(counter, v2, c2, n_spin)
!          end do
!          write(99, *) v1, c1 - val, v2, c2 - val, temp_4ks
!          coulomb_4ks_vcvc(v1, c1 - val, v2, c2 - val) = temp_4ks
!        end do
!      end do
!    end do
!  end do
!  close(99)
!  open(unit = 99, file = "coulomb_4ks_vccv")
!  do v1 = 1, val
!    do c1 = val + 1, n_states
!      do c2 = val + 1, n_states
!        do v2 = 1, val
!          temp_4ks = 0
!          do counter = 1, n_basbas
!            temp_4ks = temp_4ks + ovlp_3ks_bse(counter, v1, c1, n_spin) *
!ovlp_3ks_bse(counter, c2, v2, n_spin)
!          end do
!          write(99, *) v1, c1 - val, c2 - val, v2,  temp_4ks
!          coulomb_4ks_vccv(v1, c1 - val, c2 - val, v2) = temp_4ks
!        end do
!      end do
!    end do
!  end do
!  close(99)
  deallocate(ovlp_4ks_temp)
  deallocate(ovlp_3ks_vc)
  deallocate(ovlp_3ks_cv)
end subroutine
