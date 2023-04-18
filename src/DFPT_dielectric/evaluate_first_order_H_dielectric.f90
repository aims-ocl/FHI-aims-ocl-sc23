!****s* FHI-aims/evaluate_first_order_H_dielectric
!  NAME
!    evaluate_first_order_H_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_H_dielectric& 
         ( n_points, partition_tab, & 
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           H_times_psi, &  
           first_order_rho,first_order_potential, & 
           dVxc_drho, &
           vsigma, v2rho2, v2rhosigma, v2sigma2, & 
           gradient_rho, first_order_gradient_rho, &
           first_order_h_sparse, i_batch )

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements for phonon_gamma
!    four terms:
!   (1) <X0| r   |X0> 
!   (2) <X0| delta_V_hartree(1) |X0> 
!   (2) <X0| dVxc/drho * rho(1) |X0>   

!
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use pbc_lists
  use localorb_io, only: use_unit
  use prune_matrix_dielectric_mod

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in)    :: wave
  real*8, dimension(n_points,n_max_compute_ham)   :: wave_t,wave_t2
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, n_points), intent(in)   :: H_times_psi

  real*8, dimension(n_points), intent(in) :: first_order_rho
  real*8, dimension(n_points), intent(in) :: first_order_potential
  !-------LDA---------------
  real*8, dimension(3,n_points), intent(in) :: dVxc_drho
  !-------GGA----------------
  real*8, dimension(3,n_points), intent(in) :: v2rho2  ! = dVxc_drho
  real*8, dimension(3,n_points), intent(in) :: vsigma  
  real*8, dimension(6,n_points), intent(in) :: v2rhosigma
  real*8, dimension(6,n_points), intent(in) :: v2sigma2 
  real*8, dimension(3 ,n_spin, n_points), intent(in) :: gradient_rho
  real*8, dimension(3, n_spin, n_points), intent(in) :: first_order_gradient_rho
 

  real*8, dimension(n_compute_c,n_compute_c) :: first_order_H_dense, H_tmp
  real*8, dimension(n_hamiltonian_matrix_size), intent(inout) :: first_order_H_sparse
  !real*8, dimension(n_hamiltonian_matrix_size) :: first_order_H_sparse_check
  integer, optional :: i_batch
  !integer, dimension(n_compute_c,n_compute_c) :: online_sp
  !  INPUTS
  !  o  n_points -- number of grid points in this grid batch
  !  o  partition_tab -- values of partition function
  !  o  n_compute_c -- basis set computed at this subroutine
  !  o  i_basis_index -- i_compute to i_basis
  !  o  wave -- wave function in unit cell
  !  o  gradient_basis_wave -- gradient wave function in unit cell
  !  o  H_times_psi           -- H*psi  
  !  o  first_order_rho       -- rho_tot(1)
  !  o  first_order_potential -- V_free_hartree(1)+delta_V_hartree(1)
  !  o  dVxc_drho             -- d Vxc/ drho
  !
  !  OUTPUT
  !  first_order_H 


  integer :: i_compute,j_compute, i_basis,j_basis, i_basis_uc, j_basis_uc
  integer :: i_cell, j_cell
  integer :: i_point, i_place, i_spin, i_coord, count1, count2
  real*8  :: point_term
  real*8 :: basis_basis, gradient_basis_basis(3),  & 
            gradient_rho_gradient_basis_basis, & 
            first_order_gradient_rho_gradient_basis_basis, & 
            first_order_sigama(n_points), prefactor1(n_points), prefactor2(n_points), prefactor3(n_points), &
            grad_contract1(n_points,n_compute_c), contract(n_points,n_compute_c), contract_t(n_compute_c, n_points)  

  real*8, external :: ddot


  if(.not.use_gga.and..not.use_hartree_fock) then ! LDA case 
!------------------------------------LDA------------------------------------------------------

! Nath 261017: Restructuration of the code in order to stream memory more efficiently and use matrix multiplications (BLAS3)
  first_order_H_dense=0d0

  wave_t=transpose(wave) ! Exchange the dimensions

  do i_compute=1,n_compute_c
    do i_point=1,n_points
!-------------------local term: <mu| Vhartree^(1)+Vxc^(1) ---------------------
! We absorb all possible i_point-dependent quantities into "contract", this way we can use matrix multiplications (BLAS3) in the following
      contract_t(i_compute, i_point)=partition_tab(i_point)*wave_t(i_point,i_compute)*&
        (first_order_potential(i_point) + dVxc_drho(1,i_point)*first_order_rho(i_point))
    enddo
  enddo

! Calculate H^(1) in dense form: H_{i_comp,j_comp}=contract^T(i_comp,i_point)*wave_t(i_point,j_comp)

  call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
               1.d0,contract_t,n_compute_c,&
               wave_t(1:n_points, 1:n_compute_c), n_points, 0.d0, first_order_H_dense,n_compute_c)

     ! Make the matrix sparse
     !print *, "====start of a batch===="
     if (present(i_batch)) then
        !first_order_H_sparse_check = first_order_H_sparse
        call prune_matrix_dielectric_sparse_from_dense(first_order_H_sparse, first_order_H_dense, &
             n_compute_c, i_basis_index, i_batch)
        
        ! do i_compute=1,n_compute_c,1
        !    do j_compute=1,n_compute_c,1
        !       i_place = source_place(i_compute, j_compute)
        !       if (i_place > 0) then
        !          first_order_H_sparse(i_place) = first_order_H_sparse(i_place) + first_order_H_dense(i_compute, j_compute)
        !       end if
        !    end do
        ! end do
     ! end if
     ! if (.true.) then
     else
        do i_compute=1,n_compute_c,1
           i_basis    = i_basis_index(i_compute)
           i_basis_uc = Cbasis_to_basis(i_basis)
           i_cell     = center_to_cell(Cbasis_to_center(i_basis))

           do j_compute=1,n_compute_c,1
              j_basis    = i_basis_index(j_compute)
              j_basis_uc = Cbasis_to_basis(j_basis)


              if(j_basis_uc <= i_basis_uc) then
                 j_cell = center_to_cell(Cbasis_to_center(j_basis))

                 do i_place = &
                      index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
                      index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

                    if( column_index_hamiltonian( i_place) == j_basis_uc) then
                       first_order_H_sparse(i_place)=first_order_H_sparse(i_place)+first_order_H_dense(i_compute,j_compute)
                    endif
                 enddo !i_place
              endif
           enddo !j_compute
        enddo !i_compute

     endif
     !print *, "=====end of a batch====="
     !print *, maxval(first_order_H_sparse)!, maxval(abs(first_order_H_sparse - first_order_H_sparse_check))!, maxval(abs(online_sp(1:n_compute_c, 1:n_compute_c) - source_place(1:n_compute_c, 1:n_compute_c)))
     ! OLDER IMPLEMENTATION

     !  do i_compute=1,n_compute_c,1
     !     i_basis    = i_basis_index(i_compute)
     !     i_basis_uc = Cbasis_to_basis(i_basis)
     !     i_cell     = center_to_cell(Cbasis_to_center(i_basis))
     !
     !  do j_compute=1,n_compute_c,1
     !     j_basis    = i_basis_index(j_compute)
     !     j_basis_uc = Cbasis_to_basis(j_basis)
     !
     !     if(j_basis_uc <= i_basis_uc) then
     !        j_cell = center_to_cell(Cbasis_to_center(j_basis))
     !
     !     do i_place = &
     !        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
     !        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
     !
     !     if( column_index_hamiltonian( i_place) == j_basis_uc)then
     !
     !       do i_point = 1, n_points, 1
     !
     !           point_term = partition_tab(i_point)    &
     !                    * wave(i_compute,i_point) * wave(j_compute,i_point)
     !
     !
     !!-------------------local term: <mu| Vhartree^(1)+Vxc^(1)|nu>---------------------
     !            first_order_H_sparse(i_place) = &
     !            first_order_H_sparse(i_place) + &
     !                  point_term*first_order_potential( i_point) + &         ! (2)
     !                  point_term*dVxc_drho(1,i_point)*first_order_rho(i_point) ! (3)  
     !
     !
     !       enddo ! i_point
     !
     !     endif ! column
     !     enddo ! i_place 
     !     endif ! j_basis_uc <i_basis_uc
     !
     !
     !  enddo  ! j_compute
     !  enddo  ! i_compute 


  !if (myid .eq.0)  print*, "CHECKSPARSE", sum(first_order_H_sparse)

  else if(use_gga) then! GGA case 
!------------------------------------GGA------------------------------------------------------
! Nath 261017: Restructuration of the code in order to stream memory more efficiently and use matrix multiplications (BLAS3)
  first_order_H_dense=0d0

  i_spin = 1
  wave_t=transpose(wave) ! Exchange indices
  grad_contract1=0d0

  ! Precomputes some quantities
  ! We absorb all possible i_point-dependent quantities into different prefactors, this way we can use matrix multiplications (BLAS3) in the following
  !-------------------local term: <mu| Vhartree^(1)+Vxc^(1) ---------------------
  do i_point =1,n_points

     first_order_sigama(i_point) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
                                      gradient_rho(1:3,i_spin,i_point),1)

     prefactor1(i_point) = partition_tab(i_point)*(first_order_potential(i_point) + &       
     v2rho2(1,i_point)*first_order_rho(i_point) + &
     v2rhosigma(1,i_point)*first_order_sigama(i_point))

     prefactor2(i_point) = partition_tab(i_point)*(2.0d0*v2rhosigma(1,i_point)*first_order_rho(i_point) + &
            2.0d0*v2sigma2(1,i_point)*first_order_sigama(i_point))

     prefactor3(i_point) = partition_tab(i_point)*2.0d0*vsigma(1,i_point)

    do i_coord =1,3
      do i_compute =1,n_compute_c
        grad_contract1(i_point,i_compute) = grad_contract1(i_point,i_compute) + & 
           !----------(3.2) gradient_basis_basis term---------------
           !----------(3.2.1) gradient_rho term--------------------
         gradient_basis_wave(i_compute,i_coord,i_point)*gradient_rho(i_coord,1,i_point)*prefactor2(i_point) &
           !----------(3.2.2) first_order_gradient_rho term--------------------
        + gradient_basis_wave(i_compute,i_coord,i_point)*first_order_gradient_rho(i_coord,1,i_point)*prefactor3(i_point)
      enddo
    enddo !i_compute
  enddo !i_point

  do i_compute=1,n_compute_c
    do i_point=1,n_points
      !----------(3.1) basis_basis term-----------------------
      wave_t2(i_point,i_compute) = wave_t(i_point,i_compute)*prefactor1(i_point)
    enddo
  enddo

  ! Computes H^(1) in dense form
  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,wave_t2,n_points,&
             wave_t,n_points,0.d0,first_order_H_dense,n_compute_c)

  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,grad_contract1,n_points,&
             wave_t,n_points,0.d0,H_tmp,n_compute_c)

  first_order_H_dense=first_order_H_dense+H_tmp+transpose(H_tmp)

  !if (myid .eq.0)  print*, "CHECKDENSE", sum(first_order_H_dense)

  ! Make the matrix sparse

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)

  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

     if(j_basis_uc <= i_basis_uc) then
        j_cell = center_to_cell(Cbasis_to_center(j_basis))

        do i_place = &
           index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
           index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

           if( column_index_hamiltonian( i_place) == j_basis_uc) then

              first_order_H_sparse(i_place)=first_order_H_sparse(i_place)+first_order_H_dense(i_compute,j_compute)

           endif

        enddo !i_place
     endif
  enddo !j_compute
  enddo !i_compute

! OLDER IMPLEMENTATION

!  i_spin = 1
!!------------------------------------GGA------------------------------------------------------
!  do i_compute=1,n_compute_c,1
!     i_basis    = i_basis_index(i_compute)
!     i_basis_uc = Cbasis_to_basis(i_basis)
!     i_cell     = center_to_cell(Cbasis_to_center(i_basis))
!
!  do j_compute=1,n_compute_c,1
!     j_basis    = i_basis_index(j_compute)
!     j_basis_uc = Cbasis_to_basis(j_basis)
!
!     if(j_basis_uc <= i_basis_uc) then
!        j_cell = center_to_cell(Cbasis_to_center(j_basis))
!
!     do i_place = &
!        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
!        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
!
!     if( column_index_hamiltonian( i_place) == j_basis_uc)then
!
!       do i_point = 1, n_points, 1
!
!           basis_basis = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
!           gradient_basis_basis(1:3) = partition_tab(i_point) * ( & 
!                                       gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + &
!                                       gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) )
!  
!           gradient_rho_gradient_basis_basis(i_spin) = ddot(3, gradient_rho(1:3,i_spin,i_point), 1, & 
!                                                       gradient_basis_basis(1:3),1)
!
!           first_order_gradient_rho_gradient_basis_basis(i_spin) = & 
!                                               ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
!                                                       gradient_basis_basis(1:3),1)
!  
!           first_order_sigama(i_spin) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
!                                              gradient_rho(1:3,i_spin,i_point),1)
!
!!-------------------local term: <mu| Vhartree^(1)+Vxc^(1)|nu>---------------------
!            first_order_H_sparse(i_place) = &
!            first_order_H_sparse(i_place) + &
!                  basis_basis*first_order_potential(i_point) + &         ! (2)
!           !----------(3.1) bassis_basis term-----------------------
!           v2rho2(1,i_point)*basis_basis*first_order_rho(i_point) + &
!           v2rhosigma(1,i_point)*basis_basis*first_order_sigama(i_spin) + &
!           !----------(3.2) gradient_basis_basis term---------------
!           !----------(3.2.1) gradient_rho term--------------------
!           2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)* & 
!                                       first_order_rho(i_point) + &
!           2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)*first_order_sigama(i_spin) + &
!           !----------(3.2.2) first_order_gradient_rho term--------------------
!           2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis(i_spin)
!
!
!       enddo ! i_point
!
!     endif ! column
!     enddo ! i_place 
!     endif ! j_basis_uc <i_basis_uc
!
!
!  enddo  ! j_compute
!  enddo  ! i_compute 

  !stop
  !if (myid .eq.0)  print*, "CHECKSPARSE", sum(first_order_H_sparse)

  else if(use_hartree_fock) then   
!------------------------------------HF------------------------------------------------------
  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)

     if(j_basis_uc <= i_basis_uc) then
        j_cell = center_to_cell(Cbasis_to_center(j_basis))

     do i_place = &
        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

     if( column_index_hamiltonian( i_place) == j_basis_uc)then

       do i_point = 1, n_points, 1

           point_term = partition_tab(i_point)    &
                    * wave(i_compute,i_point) * wave(j_compute,i_point)


!-------------------local term: <mu| Vhartree^(1)+Vxc^(1)|nu>---------------------
            first_order_H_sparse(i_place) = &
            first_order_H_sparse(i_place) + &
                  point_term*first_order_potential( i_point)          ! (2)


       enddo ! i_point

     endif ! column
     enddo ! i_place 
     endif ! j_basis_uc <i_basis_uc


  enddo  ! j_compute
  enddo  ! i_compute 

  else 

  
    write(use_unit,'(2X,A)') "DFPT_dielectric is only support LDA, GGA, HF at present." 
    stop

  end if ! LDA, GGA, HF 
end subroutine evaluate_first_order_H_dielectric
!******
