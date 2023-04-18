!****s* FHI-aims/evaluate_first_order_H_polarizability
!  NAME
!    evaluate_first_order_H_polarizability
!  SYNOPSIS

subroutine evaluate_first_order_H_polarizability & 
          (first_order_H,n_points, &
           partition_tab, grid_coord,       &
           H_times_psi, n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           first_order_rho,v_hartree_gradient,dVxc_drho,& 
           vsigma, v2rho2, v2rhosigma, v2sigma2, & 
           gradient_rho, first_order_gradient_rho )

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements.
!    five terms:
!   (1) <X0| -r                 |X0>   ----------
!   (2) <X0| Int{rho(1)/|r-r'|} |X0>     Vscf(1)
!   (3) <X0| dVxc/drho * rho(1) |X0>   ----------

!  shanghui,2012.03.05
!  shanghui,2012.04.23
!  shanghui,2012.05.29 : to complie 
!  shanghui, 2013.12.09 for polarizability
 
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab
  real*8, dimension(3, n_points), intent(in) :: grid_coord

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points,n_spin), intent(in) :: H_times_psi

  real*8, dimension(3, n_spin, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_points), intent(in) :: v_hartree_gradient

  !-------LDA----------------
  real*8, dimension(3,n_points), intent(in) :: dVxc_drho

  !-------GGA----------------
  real*8, dimension(3,n_points), intent(in) :: v2rho2  ! = dVxc_drho
  real*8, dimension(3,n_points), intent(in) :: vsigma  
  real*8, dimension(6,n_points), intent(in) :: v2rhosigma
  real*8, dimension(6,n_points), intent(in) :: v2sigma2 
  real*8, dimension(3 ,n_spin, n_points), intent(in) :: gradient_rho
  real*8, dimension(3, 3, n_spin, n_points), intent(in) :: first_order_gradient_rho
 
  real*8, dimension(3,n_basis,n_basis,n_spin), intent(inout) :: first_order_H
!  INPUTS
!  o  rho -- electron density
!  o  wave --
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  first_order_H 


  integer :: i_point, i_coords,i_basis,j_basis,i_compute,j_compute,i_spin,i_coords2
  real*8 :: point_term

  real*8 :: basis_basis, gradient_basis_basis(3),  & 
            gradient_rho_gradient_basis_basis(n_spin), & 
            first_order_gradient_rho_gradient_basis_basis(3), & 
            first_order_sigama(3) ,&
            first_order_sigama2(n_points,n_spin), prefactor1(n_points), prefactor2(n_points), prefactor3(n_points)
            

  real*8, external :: ddot

  real*8 :: contract(n_points,n_compute_c), wave_t(n_points,n_max_compute_ham),tmp_H(n_compute_c,n_compute_c),&
            tmp_H2(n_compute_c,n_compute_c), contract_grad(n_points,n_compute_c)
            


  if(.not.use_gga.and..not.use_hartree_fock) then ! LDA case 
!------------------------------------LDA------------------------------------------------------
  if(n_spin.eq.1) then 
     i_spin = 1 

  wave_t=transpose(wave) ! Exchange dimensions
  do i_coords=1,3
    tmp_H=0.d0

    do i_compute=1,n_compute_c
      do i_point=1,n_points
          contract(i_point,i_compute)=partition_tab(i_point)*wave_t(i_point,i_compute)*&
          (-grid_coord(i_coords,i_point)+v_hartree_gradient(i_coords,i_point) &
           +dVxc_drho(i_spin,i_point)*first_order_rho(i_coords,i_spin,i_point))
      enddo
    enddo

    call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
               1.d0,contract,n_points,&
               wave_t,n_points,0.d0,tmp_H,n_compute_c)

    do j_compute=1,n_compute_c,1
       j_basis=i_basis_index(j_compute)
    do i_compute=1,n_compute_c,1
       i_basis=i_basis_index(i_compute)
  
         !do i_coords = 1, 3 
  
            first_order_H(i_coords,i_basis,j_basis,i_spin) =    &
            first_order_H(i_coords,i_basis,j_basis,i_spin) + tmp_H(i_compute,j_compute) 
  
         !enddo
  
    enddo
    enddo

  enddo !i_coords



!   print*, 'SUM0', sum(first_order_H)

! Old way of calculating first_order_H (slow)

!   do i_compute=1,n_compute_c,1
!   do j_compute=1,n_compute_c,1
!  
!      i_basis=i_basis_index(i_compute)
!      j_basis=i_basis_index(j_compute)
!  
!   do i_point = 1, n_points, 1
!  
!      point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
!  
!        do i_coords = 1, 3 
!  
!           first_order_H(i_coords,i_basis,j_basis,i_spin) =    &
!           first_order_H(i_coords,i_basis,j_basis,i_spin) +    & 
!           point_term*(-grid_coord(i_coords, i_point))   +    & !(1) 
!           point_term*v_hartree_gradient(i_coords,i_point)    +    & !(2)
!           point_term*dVxc_drho(i_spin,i_point)*first_order_rho(i_coords,i_spin,i_point)  !(3)
!  
!         enddo
!  
!   end do ! n_points
!  
!   enddo
!   enddo
!
!   print*, 'SUM1', sum(first_order_H)

  else !n_spin=2
  
   do i_compute=1,n_compute_c,1
   do j_compute=1,n_compute_c,1
  
      i_basis=i_basis_index(i_compute)
      j_basis=i_basis_index(j_compute)
  
   do i_point = 1, n_points, 1
  
      point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
  
        do i_coords = 1, 3 
 
           ! i_spin = 1
           ! up= dK2/(drho_u drho_u) drho_u/dE + dK2/(drho_u drho_d) drho_d/dE 
           first_order_H(i_coords,i_basis,j_basis,1) =    &
           first_order_H(i_coords,i_basis,j_basis,1) +    & 
           point_term*(-grid_coord(i_coords, i_point))   +    & !(1) 
           point_term*v_hartree_gradient(i_coords,i_point)    +    & !(2)
           point_term*dVxc_drho(1,i_point)*first_order_rho(i_coords,1,i_point) + & !(3)
           point_term*dVxc_drho(2,i_point)*first_order_rho(i_coords,2,i_point)

           ! i_spin = 2
           ! down= dK2/(drho_d drho_d) drho_d/dE + dK2/(drho_u drho_d) drho_u/dE 
           first_order_H(i_coords,i_basis,j_basis,2) =    &
           first_order_H(i_coords,i_basis,j_basis,2) +    & 
           point_term*(-grid_coord(i_coords, i_point))   +    & !(1) 
           point_term*v_hartree_gradient(i_coords,i_point)    +    & !(2)
           point_term*dVxc_drho(3,i_point)*first_order_rho(i_coords,2,i_point) + & !(3)
           point_term*dVxc_drho(2,i_point)*first_order_rho(i_coords,1,i_point)
  
         enddo
  
   end do ! n_points
  
   enddo
   enddo

  endif !n_spin

  else if(use_gga) then! GGA case 
!------------------------------------GGA------------------------------------------------------
 
  if(n_spin.eq.1) then 
     i_spin = 1

! 28.06.19 Nath: new structure for faster evaluations

  wave_t=transpose(wave) ! Exchange dimensions
  do i_coords = 1, 3
    tmp_H = 0.d0
    tmp_H2 = 0.d0
    contract_grad = 0.0

  do i_point = 1, n_points
    ! Nath: This use of ddot is very slow here...must depend on the position of the argument which is multiplied (here in 2nd position)
    !first_order_sigama(i_spin) = 2.0d0*ddot(3, first_order_gradient_rho(i_coords,1:3,i_spin,i_point), 1, &
    !                                   gradient_rho(1:3,i_spin,i_point),1)
    first_order_sigama2(i_point,i_spin) = 0.0
    do i_coords2=1,3
      first_order_sigama2(i_point,i_spin) = first_order_sigama2(i_point,i_spin) &
       +2.0d0*first_order_gradient_rho(i_coords,i_coords2,i_spin,i_point)*gradient_rho(i_coords2,i_spin,i_point)
    enddo

    ! Put together i_point-dependent factors
    prefactor1(i_point) = &
    partition_tab(i_point) *&
    (- grid_coord(i_coords, i_point) + v_hartree_gradient(i_coords,i_point) &
     + v2rho2(1,i_point)*first_order_rho(i_coords,i_spin, i_point) &
     + v2rhosigma(1,i_point)*first_order_sigama2(i_point,i_spin) )

    prefactor2(i_point) = partition_tab(i_point)*(2.0d0*v2rhosigma(1,i_point)*first_order_rho(i_coords,i_spin,i_point) &
                          + 2.0d0*v2sigma2(1,i_point)*first_order_sigama2(i_point,i_spin))

    prefactor3(i_point) = partition_tab(i_point)*2.0d0*vsigma(1,i_point)

    do i_compute = 1, n_compute_c

      contract(i_point,i_compute) = prefactor1(i_point)*wave(i_compute,i_point)

      !----------(3.2) gradient_basis_basis term---------------
      do i_coords2 = 1, 3
        !----------(3.2.1) gradient_rho term--------------------
        contract_grad(i_point,i_compute) = contract_grad(i_point,i_compute)  &
        + prefactor2(i_point)*gradient_basis_wave(i_compute,i_coords2,i_point)*gradient_rho(i_coords2,i_spin,i_point) &
        !----------(3.2.2) first_order_gradient_rho term--------------------
        + prefactor3(i_point)*gradient_basis_wave(i_compute,i_coords2,i_point) &
          *first_order_gradient_rho(i_coords,i_coords2,i_spin,i_point)
      enddo


    enddo ! i_compute
  enddo ! i_point

  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,contract,n_points,&
             wave_t,n_points,0.d0,tmp_H,n_compute_c)

  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,contract_grad,n_points,&
             wave_t,n_points,0.d0,tmp_H2,n_compute_c)

  tmp_H = tmp_H + tmp_H2 + transpose(tmp_H2)

  ! "Reconstruct" matrix
  do i_compute=1,n_compute_c,1
    i_basis=i_basis_index(i_compute)
    do j_compute=1,n_compute_c,1
      j_basis=i_basis_index(j_compute)

      first_order_H(i_coords,i_basis,j_basis,i_spin) =    &
      first_order_H(i_coords,i_basis,j_basis,i_spin) + tmp_H(i_compute,j_compute) 

    enddo
  enddo

enddo ! i_coords

! ENd new structure


!!!!! BEGIN OLDER METHOD (slower)
!
!! Nath: I made a few changes (order of loops and so on) to make this GGA part faster, but only did it for n_spin=1. n_spin=2 still needs to be done.
!
!   do i_point = 1, n_points, 1
!     do i_compute=1,n_compute_c,1
!       i_basis=i_basis_index(i_compute)
!       do j_compute=1,n_compute_c,1
!         j_basis=i_basis_index(j_compute)
!
!         basis_basis = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
!
!         gradient_basis_basis(1:3) = partition_tab(i_point) * ( & 
!                                     gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) +&
!                                     gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) )
!           
!         gradient_rho_gradient_basis_basis(i_spin)=0.0
!         do i_coords=1,3
!           gradient_rho_gradient_basis_basis(i_spin) = gradient_rho_gradient_basis_basis(i_spin) &
!                             +gradient_rho(i_coords,i_spin,i_point)*gradient_basis_basis(i_coords)
!         enddo
!  
!         do i_coords = 1, 3 
! 
!           ! Nath: This use of ddot is very slow here...must depend on the position of the argument which is multiplied (here in 2nd position)
!           !first_order_gradient_rho_gradient_basis_basis(i_spin) = & 
!           !                                    ddot(3, first_order_gradient_rho(i_coords,1:3,i_spin,i_point), 1, &
!           !                                            gradient_basis_basis(1:3),1)
!           first_order_gradient_rho_gradient_basis_basis(i_spin) = 0.0
!           do i_coords2=1,3
!             first_order_gradient_rho_gradient_basis_basis(i_spin) = first_order_gradient_rho_gradient_basis_basis(i_spin) &
!           + first_order_gradient_rho(i_coords,i_coords2,i_spin,i_point)*gradient_basis_basis(i_coords2)
!           !  first_order_gradient_rho_gradient_basis_basis(i_spin) = first_order_gradient_rho_gradient_basis_basis(i_spin) &
!           !+ first_order_gradient_rho_2(i_point,i_coords,i_coords2,i_spin)*gradient_basis_basis(i_coords2)
!           enddo
!              
!  
!           ! Nath: This use of ddot is very slow...
!           !first_order_sigama(i_spin) = 2.0d0*ddot(3, first_order_gradient_rho(i_coords,1:3,i_spin,i_point), 1, &
!           !                                   gradient_rho(1:3,i_spin,i_point),1)
!           first_order_sigama(i_spin) = 0.0
!           do i_coords2=1,3
!             first_order_sigama(i_spin) = first_order_sigama(i_spin) &
!              +2.0d0*first_order_gradient_rho(i_coords,i_coords2,i_spin,i_point)*gradient_rho(i_coords2,i_spin,i_point)
!           enddo
!  
!           first_order_H(i_coords,i_basis,j_basis,i_spin) =    &
!           first_order_H(i_coords,i_basis,j_basis,i_spin) +    & 
!           basis_basis*(-grid_coord(i_coords, i_point))   +    & !(1) 
!           basis_basis*v_hartree_gradient(i_coords,i_point)    +    & !(2)
!           !----------(3.1) bassis_basis term-----------------------
!           v2rho2(1,i_point)*basis_basis*first_order_rho(i_coords,i_spin, i_point) + &
!           v2rhosigma(1,i_point)*basis_basis*first_order_sigama(i_spin) + &
!           !----------(3.2) gradient_basis_basis term---------------
!           !----------(3.2.1) gradient_rho term--------------------
!           2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)* & 
!                                       first_order_rho(i_coords,i_spin,i_point) + &
!           2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)*first_order_sigama(i_spin) + &
!           !----------(3.2.2) first_order_gradient_rho term--------------------
!           2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis(i_spin)
!         
!         enddo ! i_coords
!  
!   enddo ! j_compute
!   enddo ! i_compute
!   end do ! i_point
!
!
!!!!! END OLDER METHOD


  else  !n_spin = 2 for GGA

   do i_compute=1,n_compute_c,1
   do j_compute=1,n_compute_c,1
  
      i_basis=i_basis_index(i_compute)
      j_basis=i_basis_index(j_compute)
  
   do i_point = 1, n_points, 1
  
         basis_basis = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
         gradient_basis_basis(1:3) = partition_tab(i_point) * ( & 
                                     gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + &
                                     gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) )
 
         do i_spin = 1,n_spin 
         gradient_rho_gradient_basis_basis(i_spin) = ddot(3, gradient_rho(1:3,i_spin,i_point), 1, & 
                                                          gradient_basis_basis(1:3),1)
         enddo  

        do i_coords = 1, 3 
 
           do i_spin = 1, n_spin 
              first_order_gradient_rho_gradient_basis_basis(i_spin) = & 
                                               ddot(3, first_order_gradient_rho(i_coords,1:3,i_spin,i_point), 1, &
                                                       gradient_basis_basis(1:3),1)
           enddo 
 
           first_order_sigama(1) = 2.0d0*ddot(3, first_order_gradient_rho(i_coords,1:3,1,i_point), 1, &
                                              gradient_rho(1:3,1,i_point),1)
           first_order_sigama(2) = ddot(3, first_order_gradient_rho(i_coords,1:3,1,i_point), 1, &
                                              gradient_rho(1:3,2,i_point),1) + & 
                                   ddot(3, first_order_gradient_rho(i_coords,1:3,2,i_point), 1, & 
                                              gradient_rho(1:3,1,i_point),1) 
           first_order_sigama(3) = 2.0d0*ddot(3, first_order_gradient_rho(i_coords,1:3,2,i_point), 1, &
                                              gradient_rho(1:3,2,i_point),1)
 
         !===================i_spin=1================================
           first_order_H(i_coords,i_basis,j_basis,1) =    &
           first_order_H(i_coords,i_basis,j_basis,1) +    & 
           basis_basis*(-grid_coord(i_coords, i_point))   +    & !(1) 
           basis_basis*v_hartree_gradient(i_coords,i_point)    +    & !(2)
           !----------(3.1) bassis_basis term-----------------------
           v2rho2(1,i_point)*basis_basis*first_order_rho(i_coords,1, i_point) + &
           v2rho2(2,i_point)*basis_basis*first_order_rho(i_coords,2, i_point) + &
           v2rhosigma(1,i_point)*basis_basis*first_order_sigama(1) + &
           v2rhosigma(2,i_point)*basis_basis*first_order_sigama(2) + &
           v2rhosigma(3,i_point)*basis_basis*first_order_sigama(3) + &
           !----------(3.2) gradient_basis_basis term---------------
           !----------(3.2.1) gradient_rho term--------------------
          ( 2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(i_coords,1,i_point) + &
                  v2rhosigma(2,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(i_coords,1,i_point)) + &
          ( 2.0d0*v2rhosigma(4,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(i_coords,2,i_point) + &
                  v2rhosigma(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(i_coords,2,i_point)) + &
          ( 2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(1) + & 
                  v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(1)) + &
          ( 2.0d0*v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(2) + &   
                  v2sigma2(4,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(2)) + &
          ( 2.0d0*v2sigma2(3,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(3) + &
                  v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(3)) + &
          !----------(3.2.2) first_order_gradient_rho term--------------------
          ( 2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis(1) + &
                  vsigma(2,i_point)*first_order_gradient_rho_gradient_basis_basis(2)) 

         !===================i_spin=2================================
           first_order_H(i_coords,i_basis,j_basis,2) =    &
           first_order_H(i_coords,i_basis,j_basis,2) +    & 
           basis_basis*(-grid_coord(i_coords, i_point))   +    & !(1) 
           basis_basis*v_hartree_gradient(i_coords,i_point)    +    & !(2)
           !----------(3.1) bassis_basis term-----------------------
           v2rho2(2,i_point)*basis_basis*first_order_rho(i_coords,1, i_point) + &
           v2rho2(3,i_point)*basis_basis*first_order_rho(i_coords,2, i_point) + &
           v2rhosigma(4,i_point)*basis_basis*first_order_sigama(1) + &
           v2rhosigma(5,i_point)*basis_basis*first_order_sigama(2) + &
           v2rhosigma(6,i_point)*basis_basis*first_order_sigama(3) + &
          !----------(3.2) gradient_basis_basis term---------------
          !----------(3.2.1) gradient_rho term--------------------
          ( 2.0d0*v2rhosigma(3,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(i_coords,1,i_point) + &
                  v2rhosigma(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(i_coords,1,i_point)) + &
          ( 2.0d0*v2rhosigma(6,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(i_coords,2,i_point) + &
                  v2rhosigma(5,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(i_coords,2,i_point)) + &
          ( 2.0d0*v2sigma2(3,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(1) + &       
                  v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(1)) + &
          ( 2.0d0*v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(2) + &       
                  v2sigma2(4,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(2)) + &
          ( 2.0d0*v2sigma2(6,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(3) + &  
                  v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(3)) + &
          !----------(3.2.2) first_order_gradient_rho term--------------------
          ( 2.0d0*vsigma(3,i_point)*first_order_gradient_rho_gradient_basis_basis(2) + &
                  vsigma(2,i_point)*first_order_gradient_rho_gradient_basis_basis(1)) 

         
         enddo
  
   end do ! n_points
  
   enddo
   enddo

  endif !n_spin 

  else if(use_hartree_fock) then   
!------------------------------------HF------------------------------------------------------
!              here I only add T + V_coulomb part, the Fock matrix is added later.
!              At present, I mean pure HF here. 
  do i_spin = 1 , n_spin
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

     do i_point = 1, n_points, 1

     point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 

       do i_coords = 1, 3 

          first_order_H(i_coords,i_basis,j_basis,i_spin) =    &
          first_order_H(i_coords,i_basis,j_basis,i_spin) +    & 
          point_term*(-grid_coord(i_coords, i_point))   +    & !(1) 
          point_term*v_hartree_gradient(i_coords,i_point)   !(2)
           
        enddo

      end do ! n_points

  enddo
  enddo
  enddo !i_spin

  else 

  
    write(use_unit,'(2X,A)') "DFPT_polarizability is only support LDA, GGA, HF at present." 
    stop

  end if ! LDA, GGA, HF 

 
end subroutine evaluate_first_order_H_polarizability
!******
