subroutine screened_coulomb_4ks(val, cond, screened_coulomb_bse, ovlp_3ks_bse, screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv)
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  implicit none
! arguments
  integer :: val, cond
  real*8, dimension(n_basbas, n_loc_prodbas) :: screened_coulomb_bse
  real*8, dimension(n_basbas,n_states,n_states,n_spin) :: ovlp_3ks_bse
  real*8, dimension(val, val, cond, cond) :: screened_coulomb_4ks_vvcc
  real*8, dimension(val, cond, cond, val) :: screened_coulomb_4ks_vccv
! local variables
  real*8 :: temp_4ks
  real*8, dimension(:, :), allocatable :: ovlp_3ks_vv, ovlp_3ks_cc
  real*8, dimension(:, :), allocatable :: ovlp_3ks_vc, ovlp_3ks_cv
  real*8, dimension(:, :), allocatable :: sc_ovlp
  real*8, dimension(:, :), allocatable :: sc_ovlp_all_states, ovlp_3ks_all_states
  real*8, dimension(:, :), allocatable :: ovlp_4ks_temp
  integer :: v1, v2, c1, c2, counter, col, row, prod_vc, prod_vv, prod_cc
! debug
  real*8 :: w_int_4ks

  prod_vv = val * val
  prod_cc = cond * cond
  prod_vc = val * cond
!  print*, ovlp_3ks_bse
!*****************sc_ovlp_all_states*****************
  allocate(sc_ovlp_all_states(n_basbas, n_states * n_states))
  sc_ovlp_all_states = 0
  allocate(ovlp_3ks_all_states(n_basbas, n_states * n_states))
  ovlp_3ks_all_states = 0
  do row = 1, n_basbas
    do c1 = 1, n_states
      do c2 = 1, n_states
        ovlp_3ks_all_states(row, (c1 - 1) * n_states + c2) = ovlp_3ks_bse(row, c1, c2, n_spin)
      end do
    end do
  end do
  call dgemm('N', 'N', n_basbas, n_states * n_states, n_basbas, 1.0d0,&
             screened_coulomb_bse, n_basbas, ovlp_3ks_all_states, n_basbas, 0.0d0, sc_ovlp_all_states, n_basbas)
!  open(unit = 99, file = 'sc_ovlp')
!!  do row = 1, n_basbas
!    row = 1
!    do c1 = 1, n_states
!      do c2 = 1, n_states
!        write(99, *) sc_ovlp_all_states(row, (c1 - 1) * n_states + c2)
!!        write(99, *) sc_ovlp_all_states(row, (c1 - val - 1) * cond + c2 - val)
!!        write(99, *) row, c1, c2, sc_ovlp(row, (c1 - val - 1) * cond + c2 - val)
!      end do
!    end do
!!  end do
!  close(99)
!***************************************************
!
  allocate(ovlp_3ks_cc(n_basbas, prod_cc))
  do row = 1, n_basbas
    do c1 = val + 1, n_states
      do c2 = val + 1, n_states
        ovlp_3ks_cc(row, (c1 - val - 1) * cond + c2 - val)  = ovlp_3ks_bse(row, c1, c2, n_spin)
      end do
    end do
  end do
!  print*, 'end initiating ovlp_3ks_vv, ovlp_3ks_cc'
!  allocate(sc_ovlp(n_basbas, n_states * n_states))
  allocate(sc_ovlp(n_basbas, prod_cc))
  sc_ovlp = 0
  call dgemm('N', 'N', n_basbas, prod_cc, n_basbas, 1.0d0,&
             screened_coulomb_bse, n_basbas, ovlp_3ks_cc, n_basbas, 0.0d0, sc_ovlp, n_basbas)
  if(allocated(ovlp_3ks_cc)) deallocate(ovlp_3ks_cc)
  print*, 'end calculating sc_ovlp'
! print to check
!  open(unit = 99, file = 'sc_ovlp')
  do row = 1, n_basbas
    do c1 = val + 1, n_states
      do c2 = val + 1, n_states
!        write(99, *) sc_ovlp(row, (c1 - val - 1) * cond + c2 - val)
!        write(99, *) row, c1, c2, sc_ovlp(row, (c1 - val - 1) * cond + c2 - val)
      end do
    end do
  end do
!  close(99)

  allocate(ovlp_3ks_vv(n_basbas, prod_vv))
  do row = 1, n_basbas
    do v1 = 1, val
      do v2 = 1, val
        ovlp_3ks_vv(row, (v1 - 1) * val + v2) = ovlp_3ks_bse(row, v1, v2, n_spin)
      end do
    end do
  end do
  allocate(ovlp_4ks_temp(prod_vv, prod_cc))
  ovlp_4ks_temp = 0
  call dgemm('T', 'N', prod_vv, prod_cc, n_basbas, 1.0d0,&
             ovlp_3ks_vv, n_basbas, sc_ovlp, n_basbas, 0.0d0, ovlp_4ks_temp, prod_vv)
  if(allocated(sc_ovlp)) deallocate(sc_ovlp)
  if(allocated(ovlp_3ks_vv)) deallocate(ovlp_3ks_vv)
  print*, 'end calculating ovlp_4ks_temp'

!  open(unit = 99, file = "screened_coulomb_4ks_vvcc_dgemm")
  screened_coulomb_4ks_vvcc = 0
  do v1 = 1, val
    do v2 = 1, val
      do c1 = 1, cond
        do c2 = 1, cond
          screened_coulomb_4ks_vvcc(v1, v2, c1, c2) = ovlp_4ks_temp((v1 - 1) * val + v2, (c1 - 1) * cond + c2) 
!          if(abs(screened_coulomb_4ks_vvcc(v1, v2, c1, c2)) > 1e-9) write(99, *) v1, v2, c1, c2, screened_coulomb_4ks_vvcc(v1, v2, c1, c2)
        end do
      end do
    end do
  end do
!  close(99)
  if(allocated(ovlp_4ks_temp)) deallocate(ovlp_4ks_temp)
  print*, 'end calculating screened_coulomb_4ks_vvcc '

!*********************************************************
  allocate(ovlp_3ks_vc(n_basbas, prod_vc))
  do row = 1, n_basbas
    do v1 = 1, val
      do c1 = val + 1, n_states
        ovlp_3ks_vc(row, (v1 - 1) * cond + c1 - val)  = ovlp_3ks_bse(row, v1, c1, n_spin)
      end do
    end do
  end do
!  print*, 'end initiating ovlp_3ks_vc'
!  ovlp_3ks_vc = 0
!  ovlp_3ks_cv = 0
  allocate(ovlp_3ks_cv(n_basbas, prod_vc))
  do row = 1, n_basbas
    do c1 = val + 1, n_states
      do v1 = 1, val
        ovlp_3ks_cv(row, (c1 - val - 1) * val + v1)  = ovlp_3ks_bse(row, c1, v1, n_spin)
      end do
    end do
  end do
!  print*, 'end initiating ovlp_3ks_cv'

  allocate(sc_ovlp(n_basbas, prod_vc)) 
!  print*, 'end allocating sc_ovlp'
  sc_ovlp = 0
  call dgemm('N', 'N', n_basbas, prod_vc, n_basbas, 1.0d0,&
             screened_coulomb_bse, n_basbas, ovlp_3ks_cv, n_basbas, 0.0d0, sc_ovlp, n_basbas)
!  print*, 'end calculating sc_ovlp'

  allocate(ovlp_4ks_temp(prod_vc, prod_vc))
  ovlp_4ks_temp = 0
  call dgemm('T', 'N', prod_vc, prod_vc, n_basbas, 1.0d0,&
             ovlp_3ks_vc, n_basbas, sc_ovlp, n_basbas, 0.0d0, ovlp_4ks_temp, prod_vc)
  if(allocated(sc_ovlp)) deallocate(sc_ovlp)
  if(allocated(ovlp_3ks_cv)) deallocate(ovlp_3ks_cv)
  if(allocated(ovlp_3ks_vc)) deallocate(ovlp_3ks_vc)

  do v1 = 1, val
    do c1 = 1, cond
      do c2 = 1, cond
        do v2 = 1, val
          screened_coulomb_4ks_vccv(v1, c1, c2, v2) = ovlp_4ks_temp((v1 - 1) * cond + c1, (c2 - 1) * val + v2)
        end do
      end do
    end do
  end do
  if(allocated(ovlp_4ks_temp)) deallocate(ovlp_4ks_temp)

!            do v1 = 1, val
!              do c2 = val + 1, n_states
!                do c1 = val + 1, n_states
!                  do v2 = 1, val
!                    w_int_4ks = 0
!                    do row = 1, n_basbas
!                      do col = 1, n_basbas
!                        w_int_4ks = w_int_4ks + ovlp_3ks_bse(row, v1, c2, n_spin)&
!                                    * screened_coulomb_bse(row, col) * ovlp_3ks_bse(col, c1, v2, n_spin)
!                      end do
!                    end do
!                    screened_coulomb_4ks_vccv(v1, c2 - val, c1 - val, v2) = w_int_4ks
!                 end do
!                end do
!              end do
!            end do
end subroutine
