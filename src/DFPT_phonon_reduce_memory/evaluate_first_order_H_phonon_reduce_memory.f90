!****s* FHI-aims/evaluate_first_order_H_phonon_reduce_memory
!  NAME
!    evaluate_first_order_H_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_H_phonon_reduce_memory& 
         ( n_points, partition_tab, & 
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, kinetic_gradient_basis_wave, & 
           H_times_psi, &  
           first_order_rho,first_order_potential, & 
           dVxc_drho, &
           coords_npoints, & 
           i_q_point, j_atom, j_coord,  &
           first_order_H_sparse)

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements for phonon_gamma
!    four terms:
!   (1) <X0| free_V_hartree(1) |X0> 
!   (2) <X0| delta_V_hartree(1) |X0> 
!   (2) <X0| dVxc/drho * rho(1) |X0>   
!   (3) <X0| Hks(0)             |X1>
!   (4) <X1| Hks(0)             |X0>

!   shanghui  2013.12.30
!
!   shanghui  2015.07. change to phonon_reduce_meory
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
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: kinetic_gradient_basis_wave   !!!DSB
  real*8, dimension(n_max_compute_ham, n_points), intent(in) :: H_times_psi

  complex*16, dimension(n_points), intent(in) :: first_order_rho
  complex*16, dimension(n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in) :: dVxc_drho

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

  real*8, external :: ddot

  real*8, allocatable ::  first_order_wave(:,:)
  real*8, allocatable ::  first_order_kinetic_wave(:,:)
  real*8, allocatable ::  contract(:,:)
  real*8, allocatable ::  first_order_H_dense(:,:)
  real*8, allocatable ::  first_order_H_dense_temp(:,:)

  allocate(first_order_wave(n_compute_c, n_points))
  allocate(first_order_kinetic_wave(n_compute_c, n_points))
  allocate(contract(n_points,n_compute_c))
  allocate(first_order_H_dense(n_compute_c,n_compute_c))
  allocate(first_order_H_dense_temp(n_compute_c,n_compute_c))

!!! !------------------[1] This is the naive loop verion---------------------------
!!!do i_compute=1,n_compute_c,1
!!!   i_basis    = i_basis_index(i_compute)
!!!   i_basis_uc = Cbasis_to_basis(i_basis)
!!!   i_cell     = center_to_cell(Cbasis_to_center(i_basis))
!!! 
!!!do j_compute=1,n_compute_c,1
!!!   j_basis    = i_basis_index(j_compute)
!!!   j_basis_uc = Cbasis_to_basis(j_basis)
!!! 
!!!   if(j_basis_uc <= i_basis_uc) then
!!!      j_cell = center_to_cell(Cbasis_to_center(j_basis))
!!! 
!!!   do i_place = &
!!!      index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
!!!      index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
!!! 
!!!   if( column_index_hamiltonian( i_place) == j_basis_uc)then
!!! 
!!!     do i_point = 1, n_points, 1
!!! 
!!!         Gr(1) = ddot(3,recip_lattice_vector(1:3,1),1, coords_npoints(1:3,i_point),1)
!!!         Gr(2) = ddot(3,recip_lattice_vector(1:3,2),1, coords_npoints(1:3,i_point),1)
!!!         Gr(3) = ddot(3,recip_lattice_vector(1:3,3),1, coords_npoints(1:3,i_point),1)
!!! 
!!!         exp_iqr(1) = exp(-(0,1)*Gr(1)*k_point_list(i_q_point,1))
!!!         exp_iqr(2) = exp(-(0,1)*Gr(2)*k_point_list(i_q_point,2))
!!!         exp_iqr(3) = exp(-(0,1)*Gr(3)*k_point_list(i_q_point,3))
!!! 
!!!         point_term = partition_tab(i_point)    &
!!!                  * wave(i_compute,i_point) * wave(j_compute,i_point)
!!! 
!!! 
!!! !-------------------(1) Hellman-Feynman term---------------------
!!!          first_order_H_sparse(i_place) = &
!!!          first_order_H_sparse(i_place) + &
!!!                point_term*first_order_potential( i_point) + &     !(1) + (2)
!!!                point_term*dVxc_drho(i_point)*first_order_rho(i_point) !(3)  
!!! 
!!! !-------------------(2) Pulay term-------------------------------    !(4) + (5)
!!!         if(j_atom.eq.Cbasis_to_atom(i_basis)) then
!!!          first_order_H_sparse(i_place) = &
!!!          first_order_H_sparse(i_place) - &
!!!          partition_tab(i_point)* &                                               
!!!          gradient_basis_wave(i_compute,j_coord,i_point)*H_times_psi(j_compute,i_point)* &
!!!          k_phase(i_cell,i_q_point)* &
!!!          exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
!!!            if(flag_rel == REL_atomic_zora)then   !DSB symmetrizing T by hand
!!!              first_order_H_sparse(i_place) = &
!!!                  first_order_H_sparse(i_place) - &
!!!                  partition_tab(i_point)* wave(j_compute,i_point)  &
!!!                  * kinetic_gradient_basis_wave(i_compute,j_coord,i_point)* &
!!!                  k_phase(i_cell,i_q_point)* &
!!!                  exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
!!!            endif
!!!
!!!         endif 
!!! 
!!!         if(j_atom.eq.Cbasis_to_atom(j_basis)) then
!!!          first_order_H_sparse(i_place) = &
!!!          first_order_H_sparse(i_place) - &
!!!          partition_tab(i_point)* &                               
!!!          gradient_basis_wave(j_compute,j_coord,i_point)*H_times_psi(i_compute,i_point)* & 
!!!          k_phase(j_cell,i_q_point)* &
!!!          exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
!!!          if(flag_rel == REL_atomic_zora)then   !DSB symmetrizing T by hand
!!!             first_order_H_sparse(i_place) = &
!!!                  first_order_H_sparse(i_place) - &
!!!                  partition_tab(i_point)* wave(i_compute,i_point)  &
!!!                  * kinetic_gradient_basis_wave(j_compute,j_coord,i_point)* &
!!!          k_phase(j_cell,i_q_point)* &
!!!          exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
!!!
!!!         endif
!!!         endif
!!! 
!!!     enddo ! i_point
!!! 
!!!   endif ! column
!!!   enddo ! i_place 
!!!   endif ! j_basis_uc <i_basis_uc
!!! 
!!! 
!!!enddo  ! j_compute
!!!enddo  ! i_compute 


!---------------------------[2] This is the updated loop version (still use loop)---------------------------
!----here we only first consider i_q_point=0, further we will add more q point---- 
!------------------------- wave^{(1)}(u,r) --------------------
! 
!   first_order_wave(1:n_compute_c,1:n_points) = 0.0d0
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
!        point_term = partition_tab(i_point)    &
!                 * wave(i_compute,i_point) * wave(j_compute,i_point)
! 
!        first_order_H_sparse(i_place) = &
!        first_order_H_sparse(i_place) + &
!              point_term*first_order_potential( i_point)  &     !(1) + (2)
!              + point_term*dVxc_drho(i_point)*first_order_rho(i_point) & !(3)  - &
!              + partition_tab(i_point)* first_order_wave(i_compute,i_point)  & 
!                                      * H_times_psi(j_compute,i_point) &
!              + partition_tab(i_point)* first_order_wave(j_compute,i_point)  &                               
!                                      * H_times_psi(i_compute,i_point)
!      enddo ! i_point
! 
!      if(flag_rel == REL_atomic_zora)then ! shanghui add for atomic_zora
!      do i_point = 1, n_points, 1
!        first_order_H_sparse(i_place) = &
!        first_order_H_sparse(i_place) + &
!                partition_tab(i_point)* first_order_kinetic_wave(i_compute,i_point)  &
!                                      * wave(j_compute,i_point) &
!              + partition_tab(i_point)* first_order_kinetic_wave(j_compute,i_point)  &
!                                      * wave(i_compute,i_point)
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

!---------------------------[3] This is the fast BLAS verion--------------------
   first_order_wave(1:n_compute_c,1:n_points) = 0.0d0
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

!-------------------(1) <i| Vhartree^(1)+Vxc^(1) |j> ---------------------
  do i_point=1,n_points
      contract(i_point,1:n_compute_c)=partition_tab(i_point)*wave(1:n_compute_c,i_point)*&
        dble(first_order_potential(i_point) + dVxc_drho(i_point)*first_order_rho(i_point))
  enddo
  first_order_H_dense = 0.0d0
  call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
             1.d0,wave(1:n_compute_c,1:n_points),n_compute_c,&
             contract,n_points,0.d0,first_order_H_dense,n_compute_c)


!-------------------(2) <i^(1)| H |j> + <i| H |j^(1)>---------------------
  do i_point=1,n_points
      contract(i_point,1:n_compute_c)=partition_tab(i_point)*H_times_psi(1:n_compute_c,i_point)
  enddo
  first_order_H_dense_temp = 0.0d0
  call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
             1.d0,first_order_wave,n_compute_c, contract,n_points, & !<i^(1)| H |j> 
             0.d0,first_order_H_dense_temp,n_compute_c)
  first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp + &
                        transpose(first_order_H_dense_temp)
  
! do i_point=1,n_points
!   contract(i_point,1:n_compute_c)=partition_tab(i_point)*first_order_wave(1:n_compute_c,i_point)
! enddo
! first_order_H_dense_temp = 0.0d0
! call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
!            1.d0, H_times_psi(1:n_compute_c,1:n_points),n_compute_c,contract,n_points, & !<i| H |j^(1)> 
!            0.d0, first_order_H_dense_temp,n_compute_c)
! first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp

!-------------------(3) zora:<k_i^(1)|j> + <i|k_j^(1)>---------------------
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

!  do i_point=1,n_points
!    contract(i_point,1:n_compute_c)=partition_tab(i_point)*first_order_kinetic_wave(1:n_compute_c,i_point)
!  enddo
!  first_order_H_dense_temp = 0.0d0
!  call dgemm("N","N",n_compute_c,n_compute_c,n_points,&
!             1.d0,wave(1:n_compute_c,1:n_points),n_compute_c, contract,n_points, & !<i|k_j^(1)> 
!             0.d0,first_order_H_dense_temp,n_compute_c)
!  first_order_H_dense = first_order_H_dense +  first_order_H_dense_temp
  endif


!----------------------dense to sparse-----------------  
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
        first_order_H_sparse(i_place)=first_order_H_sparse(i_place)+ (1.0d0,0.0d0)*first_order_H_dense(i_compute,j_compute)
    
    endif ! column
    enddo ! i_place 
    endif ! j_basis_uc <i_basis_uc


 enddo  ! j_compute
 enddo  ! i_compute 



 deallocate(first_order_wave)
 deallocate(first_order_kinetic_wave)
 deallocate(contract)
 deallocate(first_order_H_dense)
 deallocate(first_order_H_dense_temp)



end subroutine evaluate_first_order_H_phonon_reduce_memory
!******
