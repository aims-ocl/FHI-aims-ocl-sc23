!****s* FHI-aims/evaluate_first_order_H_GGA_phonon_reduce_memory
!  NAME
!    evaluate_first_order_H_GGA_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_H_GGA_phonon_reduce_memory& 
         ( n_points, partition_tab, & 
           n_compute_c, i_basis_index, index_hessian, & 
           wave, gradient_basis_wave, kinetic_gradient_basis_wave, &
           hessian_basis_wave, & 
           H_times_psi, &  
           first_order_rho,first_order_potential, & 
           vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, &
           gradient_rho, first_order_gradient_rho, &
           coords_npoints, & 
           i_q_point, j_atom, j_coord, &
           first_order_H_sparse )

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements for phonon_reduce_meory
!    five terms:
!   (1) <X0| free_V_hartree(1) |X0> 
!   (2) <X0| delta_V_hartree(1) |X0> 
!   (3) <X0| dVxc/drho * rho(1) |X0>  + GGA_gradient_term 
!   (4) <X0| Hks(0)             |X1>
!   (5) <X1| Hks(0)             |X0>

!   shanghui  2017 @ Berlin
!
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use pbc_lists

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  integer , intent(in) :: index_hessian(3,3)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: kinetic_gradient_basis_wave   !!!DSB
  real*8, dimension(n_max_compute_ham, 6, n_points),intent(in) :: hessian_basis_wave
  real*8, dimension(n_max_compute_ham, n_points), intent(in) :: H_times_psi

  complex*16, dimension(n_points), intent(in) :: first_order_rho
  complex*16, dimension(n_points), intent(in) :: first_order_potential
  real*8, dimension(n_spin,n_points), intent(in) :: vrho
  real*8, dimension(3,n_points), intent(in) :: vsigma
  real*8, dimension(3,n_points), intent(in) :: v2rho2  ! = dVxc_drho
  real*8, dimension(6,n_points), intent(in) :: v2rhosigma
  real*8, dimension(6,n_points), intent(in) :: v2sigma2
  real*8, dimension(3,n_points), intent(in) :: gradient_rho
  complex*16, dimension(3,n_points), intent(in) :: first_order_gradient_rho


  real*8, dimension(3,n_points), intent(in) :: coords_npoints
  integer , intent(in) :: i_q_point
  integer , intent(in) :: j_atom
  integer , intent(in) :: j_coord


  complex*16, dimension(n_hamiltonian_matrix_size), intent(inout) :: first_order_H_sparse

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
  integer :: i_point, i_place
  real*8  :: Gr(3), point_term
  complex*16  ::  exp_iqr(3)
  integer :: i_coord

  real*8, external :: ddot

! !---------------------------[1] This is the naive loop version---------------------------
!  real*8 :: basis_basis, gradient_basis_basis(3),  &
!            gradient_rho_gradient_basis_basis, &
!            first_order_gradient_rho_gradient_basis_basis, &
!            first_order_sigama, &
!            first_order_gradient_basis_basis(3),  &
!            gradient_rho_first_order_gradient_basis_basis
! 
!  real*8, allocatable ::  first_order_wave(:,:)
!  real*8, allocatable ::  first_order_kinetic_wave(:,:)
!  real*8, allocatable ::  first_order_gradient_wave(:,:,:)
! 
!  allocate(first_order_wave(n_compute_c, n_points))
!  allocate(first_order_kinetic_wave(n_compute_c, n_points))
!  allocate(first_order_gradient_wave(n_compute_c,3, n_points))
! 
! !----here we only first consider i_q_point=0, further we will add more q point---- 
! !------------------------- wave^{(1)}(u,r) --------------------
! 
!   first_order_wave(1:n_compute_c,1:n_points) = 0.0d0
!   ! first_order_wave = delta(i_compute,j_atom) d wave(i_compute)/dR(j_atom)
!   do i_compute=1,n_compute_c
!      i_basis=i_basis_index(i_compute)
!      if(Cbasis_to_atom(i_basis).eq.j_atom) then
!        first_order_wave(i_compute,1:n_points) = -gradient_basis_wave(i_compute,j_coord,1:n_points)
!      endif
!   enddo
! 
!   first_order_kinetic_wave(1:n_compute_c,1:n_points) = 0.0d0
!   do i_compute=1,n_compute_c
!      i_basis=i_basis_index(i_compute)
!      if(Cbasis_to_atom(i_basis).eq.j_atom) then
!        first_order_kinetic_wave(i_compute,1:n_points) = -kinetic_gradient_basis_wave(i_compute,j_coord,1:n_points)
!      endif
!   enddo
! 
!   first_order_gradient_wave(1:n_compute_c,1:3,1:n_points) = 0.0d0
!   ! first_order_gradient_wave = delta(i_compute,j_atom) d gradient_wave(i_compute)/dR(j_atom)
!   do i_coord = 1, 3
!   do i_compute=1,n_compute_c
!      i_basis=i_basis_index(i_compute)
!      if(Cbasis_to_atom(i_basis).eq.j_atom) then 
!        first_order_gradient_wave(i_compute,i_coord,1:n_points) = &
!        -hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),1:n_points)
!      endif
!   enddo
!   enddo 
! 
! do i_compute=1,n_compute_c,1
!    i_basis    = i_basis_index(i_compute)
!    i_basis_uc = Cbasis_to_basis(i_basis)
!    i_cell     = center_to_cell(Cbasis_to_center(i_basis))
! 
! do j_compute=1,n_compute_c,1
!    j_basis    = i_basis_index(j_compute)
!    j_basis_uc = Cbasis_to_basis(j_basis)
! 
!    if(j_basis_uc <= i_basis_uc) then
!       j_cell = center_to_cell(Cbasis_to_center(j_basis))
! 
!    do i_place = &
!       index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
!       index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
! 
!    if( column_index_hamiltonian( i_place) == j_basis_uc)then
! 
!      do i_point = 1, n_points, 1
!         basis_basis =  wave(i_compute,i_point) * wave(j_compute,i_point)
!         gradient_basis_basis(1:3) = gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + &
!                                     gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) 
! 
!         gradient_rho_gradient_basis_basis = ddot(3, gradient_rho(1:3,i_point), 1, &
!                                                     gradient_basis_basis(1:3),1)
!         first_order_gradient_rho_gradient_basis_basis = &
!                                         ddot(3, dble(first_order_gradient_rho(1:3,i_point)), 1, &
!                                                 gradient_basis_basis(1:3),1)
!         first_order_sigama = 2.0d0*ddot(3, dble(first_order_gradient_rho(1:3,i_point)), 1, &
!                                            gradient_rho(1:3,i_point),1)
! 
!         do i_coord = 1,3 !loop over gradient coord
!         first_order_gradient_basis_basis(i_coord) = &
!                                 first_order_gradient_wave(i_compute,i_coord,i_point) &
!                                 *wave(j_compute,i_point)  &
!                                +first_order_wave(i_compute,i_point) &
!                                 *gradient_basis_wave(j_compute,i_coord,i_point) & 
!                                +first_order_gradient_wave(j_compute,i_coord,i_point) &
!                                 *wave(i_compute,i_point)  &
!                                +first_order_wave(j_compute,i_point) &
!                                 *gradient_basis_wave(i_compute,i_coord,i_point)  
!         enddo
!         gradient_rho_first_order_gradient_basis_basis = &
!                             ddot(3, gradient_rho(1:3,i_point), 1, &
!                                     first_order_gradient_basis_basis(1:3),1)
! 
! 
! 
!        first_order_H_sparse(i_place) = &
!        first_order_H_sparse(i_place) + &
!          partition_tab(i_point)*( & ! partition_tab only write once here.
!       !----------(1) + (2)------------------------------------
!          basis_basis*first_order_potential( i_point) + &     
!       !----------(3.1) bassis_basis term-----------------------
!          v2rho2(1,i_point)*basis_basis*first_order_rho(i_point)   + &
!          v2rhosigma(1,i_point)*basis_basis*first_order_sigama + &
!        !----------(3.2) gradient_basis_basis term---------------
!        !----------(3.2.1) gradient_rho term--------------------
!          2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis &
!                                     *first_order_rho(i_point) + &
!          2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis & 
!                                   *first_order_sigama + &
!        !----------(3.2.2) first_order_gradient_rho term--------------------
!          2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis + & 
!        !----------(4) + (5) -----------------------------------------
!          first_order_wave(i_compute,i_point)*H_times_psi(j_compute,i_point) + &
!          first_order_wave(i_compute,i_point)*vrho(1,i_point)*wave(j_compute,i_point) + &
!          first_order_wave(j_compute,i_point)*H_times_psi(i_compute,i_point) + &    
!          first_order_wave(j_compute,i_point)*vrho(1,i_point)*wave(i_compute,i_point) + &
!          2.0d0*vsigma(1,i_point)*gradient_rho_first_order_gradient_basis_basis  & 
!          ) 
! 
!      enddo ! i_point
! 
!      if(flag_rel == REL_atomic_zora)then ! shanghui add for atomic_zora
!      do i_point = 1, n_points, 1
!        first_order_H_sparse(i_place) = &
!        first_order_H_sparse(i_place) + &
!           partition_tab(i_point)* first_order_kinetic_wave(i_compute,i_point)  &
!                                 * wave(j_compute,i_point) &
!         + partition_tab(i_point)* first_order_kinetic_wave(j_compute,i_point)  &
!                                 * wave(i_compute,i_point)
!      enddo  
!      endif
! 
!    endif ! column
!    enddo ! i_place 
!    endif ! j_basis_uc <i_basis_uc
! 
! 
! enddo  ! j_compute
! enddo  ! i_compute 
! 
! deallocate(first_order_wave)
! deallocate(first_order_kinetic_wave)
! deallocate(first_order_gradient_wave)


!!---------------------------[2] This is the BLAS version---------------------------
  real*8 :: first_order_sigama, prefactor_1, prefactor_2
 
  real*8, allocatable ::  first_order_wave(:,:)
  real*8, allocatable ::  first_order_kinetic_wave(:,:)
  real*8, allocatable ::  first_order_gradient_wave(:,:,:)
  real*8, allocatable ::  contract(:,:)
  real*8, allocatable ::  first_order_H_dense(:,:)
  real*8, allocatable ::  first_order_H_dense_temp(:,:)

 
  allocate(first_order_wave(n_compute_c, n_points))
  allocate(first_order_kinetic_wave(n_compute_c, n_points))
  allocate(first_order_gradient_wave(n_compute_c,3, n_points))
  allocate(contract(n_points,n_compute_c))
  allocate(first_order_H_dense(n_compute_c,n_compute_c))
  allocate(first_order_H_dense_temp(n_compute_c,n_compute_c))


    first_order_wave(1:n_compute_c,1:n_points) = 0.0d0
    ! first_order_wave = delta(i_compute,j_atom) d wave(i_compute)/dR(j_atom)
    do i_compute=1,n_compute_c
       i_basis=i_basis_index(i_compute)
       if(Cbasis_to_atom(i_basis).eq.j_atom) then
         first_order_wave(i_compute,1:n_points) = -gradient_basis_wave(i_compute,j_coord,1:n_points)
       endif
    enddo
  
    first_order_kinetic_wave(1:n_compute_c,1:n_points) = 0.0d0
    do i_compute=1,n_compute_c
       i_basis=i_basis_index(i_compute)
       if(Cbasis_to_atom(i_basis).eq.j_atom) then
         first_order_kinetic_wave(i_compute,1:n_points) = -kinetic_gradient_basis_wave(i_compute,j_coord,1:n_points)
       endif
    enddo
  
    first_order_gradient_wave(1:n_compute_c,1:3,1:n_points) = 0.0d0
    ! first_order_gradient_wave = delta(i_compute,j_atom) d gradient_wave(i_compute)/dR(j_atom)
    do i_coord = 1, 3
    do i_compute=1,n_compute_c
       i_basis=i_basis_index(i_compute)
       if(Cbasis_to_atom(i_basis).eq.j_atom) then 
         first_order_gradient_wave(i_compute,i_coord,1:n_points) = &
         -hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),1:n_points)
       endif
    enddo
    enddo 


 !----------(1) + (2) + (3.1) bassis_basis term------------------------------------
   do i_point=1,n_points
      first_order_sigama = 2.0d0*ddot(3, dble(first_order_gradient_rho(1:3,i_point)), 1, &
                                            gradient_rho(1:3,i_point),1)
   
      contract(i_point,1:n_compute_c)=partition_tab(i_point)*wave(1:n_compute_c,i_point)*&
        dble( first_order_potential(i_point)  & 
           + v2rho2(1,i_point)*first_order_rho(i_point)  & 
           + v2rhosigma(1,i_point)*first_order_sigama )
    enddo
    first_order_H_dense = 0.0d0
    call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
               1.d0,wave(1:n_compute_c,1:n_points),n_compute_c,&
               contract,n_points,0.d0,first_order_H_dense,n_compute_c)


  
 !----------(3.2) gradient_basis_basis term------------------------------------
   contract = 0.0d0
   do i_point =1,n_points
      first_order_sigama = 2.0d0*ddot(3, dble(first_order_gradient_rho(1:3,i_point)), 1, &
                                             gradient_rho(1:3,i_point),1)
      prefactor_1 = partition_tab(i_point)*( &
          2.0d0*v2rhosigma(1,i_point)*dble(first_order_rho(i_point)) + &
          2.0d0*v2sigma2(1,i_point)*first_order_sigama)
      prefactor_2 = partition_tab(i_point)* & 
          2.0d0*vsigma(1,i_point)
   
      do i_coord =1,3
        do i_compute =1,n_compute_c
           contract(i_point,i_compute) = contract(i_point,i_compute) + &
        !----------(3.2.1) gradient_rho term--------------------
              prefactor_1 & 
             *gradient_rho(i_coord,i_point) & 
             *gradient_basis_wave(i_compute,i_coord,i_point) &
        !----------(3.2.2) first_order_gradient_rho term--------------------
             +prefactor_2  &
             *first_order_gradient_rho(i_coord,i_point) &
             *gradient_basis_wave(i_compute,i_coord,i_point)
        enddo !i_compute
      enddo !i_coord

    enddo !i_point
   
    first_order_H_dense_temp = 0.0d0
    call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
               1.d0,wave(1:n_compute_c,1:n_points),n_compute_c,&
               contract,n_points,0.d0,first_order_H_dense_temp,n_compute_c)

    first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp + &
                        transpose(first_order_H_dense_temp)

 
 !----------(4)+(5) Pulay term------------------------------------
    contract = 0.0d0
    do i_point=1,n_points
       contract(i_point,1:n_compute_c) = partition_tab(i_point)*( & 
                                         H_times_psi(1:n_compute_c,i_point) &
                                       + vrho(1,i_point)*wave(1:n_compute_c,i_point) ) 
       prefactor_1 = partition_tab(i_point)*2.0d0*vsigma(1,i_point) 

       do i_coord =1,3
       do i_compute =1,n_compute_c
           contract(i_point,i_compute) = contract(i_point,i_compute) + &
        !----------(4) + (5) first_order_wave--------------------
              prefactor_1 & 
             *gradient_rho(i_coord,i_point) & 
             *gradient_basis_wave(i_compute,i_coord,i_point) 
       enddo !i_compute
       enddo !i_coord

    enddo

    first_order_H_dense_temp = 0.0d0
    call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
               1.d0,first_order_wave,n_compute_c, contract,n_points, &  
               0.d0,first_order_H_dense_temp,n_compute_c)
    first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp + &
                          transpose(first_order_H_dense_temp)
  
 
    contract = 0.0d0
    do i_point=1,n_points
       prefactor_1 = partition_tab(i_point)*2.0d0*vsigma(1,i_point) 

       do i_coord =1,3
       do i_compute =1,n_compute_c
           contract(i_point,i_compute) = contract(i_point,i_compute) + &
        !----------(4) + (5) wave--------------------
              prefactor_1 & 
             *gradient_rho(i_coord,i_point) & 
             *first_order_gradient_wave(i_compute,i_coord,i_point) 
       enddo !i_compute
       enddo !i_coord

    enddo

    first_order_H_dense_temp = 0.0d0
    call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
               1.d0,wave(1:n_compute_c,1:n_points),n_compute_c, contract,n_points, &  
               0.d0,first_order_H_dense_temp,n_compute_c)
    first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp + &
                          transpose(first_order_H_dense_temp)
   

  !-------------------(6) zora:<k_i^(1)|j> + <i|k_j^(1)>---------------------
  if(flag_rel == REL_atomic_zora)then ! shanghui add for atomic_zora
  do i_point=1,n_points
    contract(i_point,1:n_compute_c)=partition_tab(i_point)*wave(1:n_compute_c,i_point)
  enddo
  first_order_H_dense_temp = 0.0d0
  call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
             1.d0,first_order_kinetic_wave ,n_compute_c,contract,n_points, & !<k_i^(1)|j> 
             0.d0,first_order_H_dense_temp,n_compute_c)
  first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp + &
                        transpose(first_order_H_dense_temp)
  endif


  !--------------dense to sparse--------------
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

              first_order_H_sparse(i_place)=first_order_H_sparse(i_place)+ & 
                       (1.0d0,0.0d0)*first_order_H_dense(i_compute,j_compute)

           endif

        enddo !i_place
     endif
  enddo !j_compute
  enddo !i_compute

 deallocate(first_order_wave)
 deallocate(first_order_kinetic_wave)
 deallocate(first_order_gradient_wave)
 deallocate(contract)
 deallocate(first_order_H_dense)
 deallocate(first_order_H_dense_temp)


end subroutine evaluate_first_order_H_GGA_phonon_reduce_memory
!******
