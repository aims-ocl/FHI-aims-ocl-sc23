  Subroutine CC_RLE()

! This subroutine evaluates the best coefficient vector using RLE algorithm

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  ! Calculate the new trial function t(m)
  Call CC_RLE_new_t()

  if (.not.(CC_RLE_flag_singular)) then
    ! Adding new trial vector t(m) to the basis set t(m)
    Call CC_RLE_adding_t()
  end if

  if ((.not.(CC_RLE_flag_singular)).and.(.not.(CC_RLE_flag_sv_overflow))) then
    ! Solve RLE equation alpha + R * tau = 0
    Call CC_RLE_solution()
  else
    ! Use newest trial funciton as the best trial function (Jacob method)
    Call CC_Jacob_sv_t()
  end if

  End Subroutine CC_RLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_RLE_new_t()

! This subroutine generates the new trial function for RLE algorithm

  Use dimensions
  Use physics
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_state,i_run
  Integer :: i_spin,iii,jjj,aaa,bbb,ainx,binx
  Double precision :: D_ia,D_ijab,int_rlt,et1,et2

!  print*,'new_t'

  ! For single excitation
  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,i_spin) - 1

    Call CC_decode_index_s(i_state,i_spin,iii,aaa)

    do i_run = 1,CC_mem_w(myid+1,i_spin)
      Call CC_inc_index_s(i_spin,iii,aaa)
!      print*,'iii,aaa',iii,aaa
      ! Calculate the element of D matrix
      Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,D_ia)

      ! t(new) = - D**-1[A+Delta*t], for single excitation, A = 0
      if (i_spin.eq.1) then
        CC_w_1a(i_run) = - CC_w_1a(i_run) / D_ia
      else
        CC_w_1b(i_run) = - CC_w_1b(i_run) / D_ia
      end if
!      print*,'i_run',i_run
!      print*,'CC_w_1a',CC_w_1a(i_run)
    end do
  end do

  ! For double excitation
  ! alpha-beta excitation
  i_state = CC_index_w(myid+1,3) - 1
  Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

  do i_run = 1, CC_mem_w(myid+1,3)
    Call CC_inc_index_ab(iii,jjj,aaa,bbb)

!    print*,'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
    ! Calculate element of D matrix
    Call CC_Calc_eigenvalue_diff(iii,aaa,1,1,et1)
    Call CC_Calc_eigenvalue_diff(jjj,bbb,2,2,et2)
    D_ijab = et1 + et2

    ! Calculate element of A vector
    Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,1,2,1,2)

    CC_w_ab(i_run) = - (int_rlt + CC_w_ab(i_run)) / D_ijab
!    print*,'i_run',i_run
!    print*,'CC_w_ab',CC_w_ab(i_run)
  end do

  ! alpha-alpha and/or beta-beta excitation
  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,3+i_spin) - 1

    do i_run = 1, CC_mem_w(myid+1,3+i_spin)

      i_state = i_state + 1
      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

      ! Calculate element of D matrix
      Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,et1)
      Call CC_Calc_eigenvalue_diff(jjj,bbb,i_spin,i_spin,et2)
      D_ijab = et1 + et2
  
      ! Calculate element of A vector
      Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)

      if (i_spin.eq.1) then  
        CC_w_2a(i_run) = - (int_rlt + CC_w_2a(i_run)) / D_ijab
      else
        CC_w_2b(i_run) = - (int_rlt + CC_w_2b(i_run)) / D_ijab
      end if
    end do
  end do

  End Subroutine CC_RLE_new_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_RLE_adding_t()

! This subroutine adds the new trial function to the basis set

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_run,i_local
  Integer :: i_spin,n_RLE_pro,i_pro
  Double precision :: vec_dot,t_norm,t_diff

  ! Calculate the difference between new trial function and the last basis function
  t_diff = 0.0D0
  do i_spin = 1, n_spin
    do i_local = 1, CC_mem_w(myid+1,i_spin)
      if (i_spin.eq.1) then
        t_diff = t_diff + (CC_w_1a(i_local) - CC_t_1a_sv(CC_RLE_ndc,i_local)) ** 2
      else
        t_diff = t_diff + (CC_w_1b(i_local) - CC_t_1b_sv(CC_RLE_ndc,i_local)) ** 2
      end if
    end do

    do i_local = 1, CC_mem_w(myid+1,3+i_spin)
      if (i_spin.eq.1) then
        t_diff = t_diff + (CC_w_2a(i_local) - CC_t_2a_sv(CC_RLE_ndc,i_local)) ** 2
      else
        t_diff = t_diff + (CC_w_2b(i_local) - CC_t_2b_sv(CC_RLE_ndc,i_local)) ** 2
      end if
    end do
  end do

  if (n_spin.eq.1) then
    t_diff = t_diff * 2.0D0
  end if

  do i_local = 1, CC_mem_w(myid+1,3)
    t_diff = t_diff + (CC_w_ab(i_local) - CC_t_ab_sv(CC_RLE_ndc,i_local)) ** 2
  end do

  Call sync_real_number(t_diff)

  t_diff = sqrt(t_diff)

  ! Linear dependence condition
  if (t_diff.lt.CC_RLE_bas_thresh) then
    CC_RLE_flag_singular = .true.
  end if

  ! increase the number of basis functions, if the new basis function is linear
  ! independent of the present basis set and the storage space is enough.
  if (.not.(CC_RLE_flag_singular)) then
    CC_RLE_ndc = CC_RLE_ndc + 1
    if (CC_RLE_ndc.gt.CC_n_RLE_sv) then
      CC_RLE_ndc = CC_n_RLE_sv
      CC_RLE_flag_sv_overflow = .true.
    end if
  end if

!  print*,'t_diff',t_diff
!  print*,'CC_RLE_ndc',CC_RLE_ndc

  ! Add trial function to basis set
  CC_t_1a_sv(CC_RLE_ndc,:) = CC_w_1a
  CC_t_2a_sv(CC_RLE_ndc,:) = CC_w_2a
  CC_t_ab_sv(CC_RLE_ndc,:) = CC_w_ab

  if (n_spin.ne.1) then
    CC_t_1b_sv(CC_RLE_ndc,:) = CC_w_1b
    CC_t_2b_sv(CC_RLE_ndc,:) = CC_w_2b
  end if

  End Subroutine CC_RLE_adding_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_RLE_solution()

! This subroutine solves the RLE equation

  Use dimensions
  Use mpi_tasks
  Use physics
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_run
  Integer :: iii,jjj,aaa,bbb,iinx,jinx,ainx,binx,kkk,lll,i_spin,errnum
  Integer :: i_t,NNN,NRHS,LDA,LDB,INFO,i_bas
  Integer , dimension(:) , allocatable :: IPIV
  Double precision :: D_element,Rkl,ak,int_rlt,t_new,et1,et2,t_best
  Double precision , dimension(:,:) , allocatable :: Rmat,tau
  Double precision , dimension(:) , allocatable :: avec

  NNN = CC_RLE_ndc - 1

  ! Form R matrix and a vector
  Allocate(Rmat(NNN,NNN),stat=errnum)
  Call check_allocation(errnum, 'Rmat in CC')

  Allocate(avec(NNN),stat=errnum)
  Call check_allocation(errnum, 'avec in CC')

  do kkk = 1, NNN
    do lll = 1, NNN
      ! Calculate Matrix element
      Rkl = 0.0D0
      ak = 0.0D0

      do i_spin = 1, n_spin
        i_state = CC_index_w(myid+1,i_spin)

        Call CC_decode_index_s(i_state,i_spin,iii,aaa)
        
        aaa = aaa - 1
        do i_run = 1, CC_mem_w(myid+1,i_spin)
          Call CC_inc_index_s(i_spin,iii,aaa)

          Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,D_element)

          Call CC_Calc_t_best(i_spin,lll,i_run,t_best,D_element,0.0D0)

          if (i_spin.eq.1) then
            Rkl = Rkl + CC_t_1a_sv(kkk,i_run) * D_element &
                      * (t_best - CC_t_1a_sv(lll+1,i_run))
          else
            Rkl = Rkl + CC_t_1b_sv(kkk,i_run) * D_element &
                      * (t_best - CC_t_1b_sv(lll+1,i_run))
          end if
          
        end do

        i_state = CC_index_w(myid+1,3+i_spin) - 1

        do i_run = 1, CC_mem_w(myid+1,3+i_spin)

          i_state = i_state + 1
          Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

          Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)

          Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,et1)
          Call CC_Calc_eigenvalue_diff(jjj,bbb,i_spin,i_spin,et2)
          D_element = et1 + et2

          Call CC_Calc_t_best(3+i_spin,lll,i_run,t_best,D_element,int_rlt)

          if (i_spin.eq.1) then
            Rkl = Rkl + CC_t_2a_sv(kkk,i_run) * D_element &
                      * (t_best - CC_t_2a_sv(lll+1,i_run)) &
                      - CC_t_2a_sv(kkk,i_run) * int_rlt

            ak = ak + CC_t_2a_sv(kkk,i_run) * int_rlt
          else
            Rkl = Rkl + CC_t_2b_sv(kkk,i_run) * D_element &
                      * (t_best - CC_t_2b_sv(lll+1,i_run)) &
                      - CC_t_2b_sv(kkk,i_run) * int_rlt

            ak = ak + CC_t_2b_sv(kkk,i_run) * int_rlt
          end if

        end do
      end do

      if (n_spin.eq.1) then
        Rkl = Rkl * 2.0D0
        ak = ak * 2.0D0
      end if

      i_state = CC_index_w(myid+1,3) - 1
      Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

      do i_run = 1, CC_mem_w(myid+1,3)
        Call CC_inc_index_ab(iii,jjj,aaa,bbb)

        Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,1,2,1,2)

        Call CC_Calc_eigenvalue_diff(iii,aaa,1,1,et1)
        Call CC_Calc_eigenvalue_diff(jjj,bbb,2,2,et2)

        D_element = et1 + et2

        Call CC_Calc_t_best(3,lll,i_run,t_best,D_element,int_rlt)

        Rkl = Rkl + CC_t_ab_sv(kkk,i_run) * D_element &
                  * (t_best - CC_t_ab_sv(lll+1,i_run)) &
                  - CC_t_ab_sv(kkk,i_run) * int_rlt

        ak = ak + CC_t_ab_sv(kkk,i_run) * int_rlt
      end do

      Rmat(kkk,lll) = Rkl
      avec(kkk) = ak
    end do
  end do

  Call sync_vector(Rmat,NNN*NNN)
  Call sync_vector(avec,NNN)

  ! Solve the linear equations a + R * tau = 0
  Allocate(tau(NNN,1),stat=errnum)
  Call check_allocation(errnum, 'tau in CC')

  do kkk = 1, NNN
    tau(kkk,1) = - avec(kkk)
  end do

  Allocate(IPIV(NNN),stat=errnum)
  Call check_allocation(errnum, 'IPIV in CC')

  Call DGESV(NNN,1,Rmat,NNN,IPIV,tau,NNN,INFO)

!  print*,'tau',tau

  ! Save the result 
  do i_run = 1, NNN
    CC_RLE_t(CC_RLE_ndc,i_run) = tau(i_run,1)
  end do

  ! The best evaluation of t vector
  Call CC_clear_t()

  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,i_spin) - 1
    Call CC_decode_index_s(i_state,i_spin,iii,aaa)

    do i_run = 1, CC_mem_w(myid+1,i_spin)
      Call CC_inc_index_s(i_spin,iii,aaa)
      ainx = aaa - CC_n_elec(i_spin)
      iinx = iii - CC_valence + 1

      Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,D_element)

      t_new = 0.0D0
      do i_bas = 1, NNN
        if (i_spin.eq.1) then
          t_new = t_new - D_element * CC_t_1a_sv(i_bas+1,i_run) * tau(i_bas,1)
        else
          t_new = t_new - D_element * CC_t_1b_sv(i_bas+1,i_run) * tau(i_bas,1)
        end if
      end do

      if (i_spin.eq.1) then
        CC_t_1a(iinx,ainx) = - t_new / D_element
      else
        CC_t_1b(iinx,ainx) = - t_new / D_element
      end if
    end do

    i_state = CC_index_w(myid+1,3+i_spin) - 1

    do i_run = 1, CC_mem_w(myid+1,3+i_spin)

      i_state = i_state + 1
      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

      Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)

      Call CC_Calc_eigenvalue_diff(iii,aaa,i_spin,i_spin,et1)
      Call CC_Calc_eigenvalue_diff(jjj,bbb,i_spin,i_spin,et2)

      D_element = et1 + et2

      t_new = 0.0D0
      do i_bas = 1, NNN
        if (i_spin.eq.1) then
          t_new = t_new - (D_element * CC_t_2a_sv(i_bas+1,i_run) + int_rlt) * tau(i_bas,1)
        else
          t_new = t_new - (D_element * CC_t_2b_sv(i_bas+1,i_run) + int_rlt) * tau(i_bas,1)
        end if
      end do

      if (i_spin.eq.1) then
        CC_t_2a(i_state) = - (t_new + int_rlt) / D_element
      else
        CC_t_2b(i_state) = - (t_new + int_rlt) / D_element
      end if

    end do    
  end do

  i_state = CC_index_w(myid+1,3) - 1
  Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

  do i_run = 1, CC_mem_w(myid+1,3)

    i_state = i_state + 1

    Call CC_inc_index_ab(iii,jjj,aaa,bbb)

    Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,1,2,1,2)

    Call CC_Calc_eigenvalue_diff(iii,aaa,1,1,et1)
    Call CC_Calc_eigenvalue_diff(jjj,bbb,2,2,et2)

    D_element = et1 + et2

    t_new = 0.0D0

    do i_bas = 1, NNN
      t_new = t_new - (D_element * CC_t_ab_sv(i_bas+1,i_run) + int_rlt) * tau(i_bas,1)
    end do

    CC_t_ab(i_state) = - (t_new + int_rlt) / D_element
  end do

  Call sync_vector(CC_t_1a,CC_n_1a)
  Call sync_vector(CC_t_ab,CC_n_ab)
  Call sync_vector(CC_t_2a,CC_n_2a)

  if (n_spin.ne.1) then
    Call sync_vector(CC_t_1b,CC_n_1b)
    Call sync_vector(CC_t_2b,CC_n_2b)
  end if

  Deallocate(IPIV)
  Deallocate(tau)
  Deallocate(avec)
  Deallocate(Rmat)

  End Subroutine CC_RLE_solution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_t_best(t_type,lll,i_local,t_best,D_element,a)
! This subroutine calculates the best evaluation of trial function for lll-th 
! iteration

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: t_type,lll,i_local,i_bas
  Double precision :: t_best,D_element,a,tmp

  if (lll.le.2) then
    if (t_type.eq.1) then
      t_best = CC_t_1a_sv(lll,i_local) 
    else if (t_type.eq.2) then
      t_best = CC_t_1b_sv(lll,i_local) 
    else if (t_type.eq.3) then
      t_best = CC_t_ab_sv(lll,i_local) 
    else if (t_type.eq.4) then
      t_best = CC_t_2a_sv(lll,i_local)
    else if (t_type.eq.5) then
      t_best = CC_t_2b_sv(lll,i_local) 
    end if
  else
    tmp = a / D_element
    t_best = - tmp
    do i_bas = 1, lll - 1
      if (t_type.eq.1) then
        t_best = t_best + CC_t_1a_sv(i_bas+1,i_local) * CC_RLE_t(lll,i_bas)
      else if (t_type.eq.2) then
        t_best = t_best + CC_t_1b_sv(i_bas+1,i_local) * CC_RLE_t(lll,i_bas)
      else if (t_type.eq.3) then
        t_best = t_best + (CC_t_ab_sv(i_bas+1,i_local) + tmp) * CC_RLE_t(lll,i_bas)
      else if (t_type.eq.4) then
        t_best = t_best + (CC_t_2a_sv(i_bas+1,i_local) + tmp) * CC_RLE_t(lll,i_bas)
      else if (t_type.eq.5) then
        t_best = t_best + (CC_t_2b_sv(i_bas+1,i_local) + tmp) * CC_RLE_t(lll,i_bas)
      end if
    end do
  end if

  End Subroutine CC_Calc_t_best

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
