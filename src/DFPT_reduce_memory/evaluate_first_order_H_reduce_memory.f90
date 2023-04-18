!****s* FHI-aims/evaluate_first_order_H
!  NAME
!    evaluate_first_order_H
!  SYNOPSIS

subroutine evaluate_first_order_H_reduce_memory( & 
           first_order_H,n_points, &
           partition_tab, dist_tab_global, dir_tab_global,       &
           n_compute_c, i_basis_index, index_hessian, & 
           wave, gradient_basis_wave, hessian_basis_wave, & 
           H_times_psi,  &
           first_order_rho,v_hartree_gradient, & 
           dVxc_drho, &
           vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, &
           gradient_rho, first_order_gradient_rho, &
           j_atom,j_coord,kinetic_gradient_basis_wave)

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements.
!    five terms:
!   (1) <X0| Z_I(R-r)/|r-R_I|^3 |X0>   ----------
!   (2) <X0| Int{rho(1)/|r-r'|} |X0>     Vscf(1)
!   (3) <X0| dVxc/drho * rho(1) |X0>   ----------
!   (4) <X0| Hks(0)             |X1>
!   (5) <X1| Hks(0)             |X0>

!  shanghui,2012.03.05
!  shanghui,2012.04.23
!  shanghui,2012.05.29 : to complie 
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
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_global
  real*8, dimension(3,n_atoms, n_points), intent(in) :: dir_tab_global

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  integer , intent(in) :: index_hessian(3,3)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, 6, n_points),intent(in) :: hessian_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi

  real*8, dimension(n_points), intent(in) :: first_order_rho
  real*8, dimension(n_points), intent(in) :: v_hartree_gradient
  !-------LDA---------------
  real*8, dimension(3,n_points), intent(in) :: dVxc_drho
  !-------GGA----------------
  real*8, dimension(n_spin,n_points), intent(in) :: vrho  
  real*8, dimension(3,n_points), intent(in) :: vsigma  
  real*8, dimension(3,n_points), intent(in) :: v2rho2  ! = dVxc_drho
  real*8, dimension(6,n_points), intent(in) :: v2rhosigma
  real*8, dimension(6,n_points), intent(in) :: v2sigma2 
  real*8, dimension(3 ,n_spin, n_points), intent(in) :: gradient_rho
  real*8, dimension(3, n_spin, n_points), intent(in) :: first_order_gradient_rho

  real*8, dimension(n_basis,n_basis), intent(inout) :: first_order_H

  integer, intent(in) :: j_atom, j_coord
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: kinetic_gradient_basis_wave   !!!DSB

!  INPUTS
!  o  rho -- electron density
!  o  wave --
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_global -- (distance to atoms)**1
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  first_order_H 


  integer :: i_point, i_basis,j_basis,i_compute,j_compute
  integer :: i_coord,i_spin
  real*8 :: point_term
  real*8 :: atomic_term
  real*8 :: basis_basis, gradient_basis_basis(3),  & 
            gradient_rho_gradient_basis_basis(n_spin), & 
            first_order_gradient_rho_gradient_basis_basis(3), & 
            first_order_sigama(3), &
            first_order_gradient_basis_basis(3),  &
            gradient_rho_first_order_gradient_basis_basis(n_spin) 

  real*8, external :: ddot

  if(.not.use_gga.and..not.use_hartree_fock) then ! LDA case 
!------------------------------------LDA------------------------------------------------------
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

  do i_point = 1, n_points, 1
     point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 

     atomic_term = -species_z(species(j_atom)) * point_term &
                      / (dist_tab_global(j_atom, i_point))**3

       first_order_H(i_basis,j_basis) =    &
       first_order_H(i_basis,j_basis) +    & 
          atomic_term * dir_tab_global(j_coord, j_atom, i_point)   +    & !(1) 
          point_term*v_hartree_gradient(i_point)   +    & !(2)
          point_term*dVxc_drho(1,i_point)*first_order_rho( i_point) !(3)

  end do ! n_points

  enddo
  enddo
  
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

       if(j_atom.eq.basis_atom(i_basis)) then 
       do i_point = 1, n_points, 1
       first_order_H(i_basis,j_basis) =    &
       first_order_H(i_basis,j_basis) -    &
           partition_tab(i_point)* &  
           gradient_basis_wave(i_compute,j_coord,i_point)*H_times_psi(j_compute,i_point) 

         if(flag_rel == REL_atomic_zora)then   
           first_order_H(i_basis,j_basis) =    &
           first_order_H(i_basis,j_basis) -    &
             partition_tab(i_point)* &
             wave(j_compute,i_point)*kinetic_gradient_basis_wave(i_compute,j_coord,i_point)
         endif

       enddo 
       endif 

       if(j_atom.eq.basis_atom(j_basis)) then 
       do i_point = 1, n_points, 1
       first_order_H(i_basis,j_basis) =    &
       first_order_H(i_basis,j_basis) -    &
           partition_tab(i_point)* &  
           gradient_basis_wave(j_compute,j_coord,i_point)*H_times_psi(i_compute,i_point) 

         if(flag_rel == REL_atomic_zora)then   
           first_order_H(i_basis,j_basis) =    &
           first_order_H(i_basis,j_basis) -    &
             partition_tab(i_point)*&
             wave(i_compute,i_point)*kinetic_gradient_basis_wave(j_compute,j_coord,i_point)
         endif
       enddo 
       endif
  enddo
  enddo

  else if(use_gga) then! GGA case 
  i_spin = 1
!------------------------------------GGA------------------------------------------------------
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

  do i_point = 1, n_points, 1
     point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 

     atomic_term = -species_z(species(j_atom)) * point_term &
                      / (dist_tab_global(j_atom, i_point))**3

     basis_basis = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
     gradient_basis_basis(1:3) = partition_tab(i_point) * ( & 
                                 gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + &
                                 gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) )
  
     gradient_rho_gradient_basis_basis(i_spin) = ddot(3, gradient_rho(1:3,i_spin,i_point), 1, & 
                                                 gradient_basis_basis(1:3),1)

     first_order_gradient_rho_gradient_basis_basis(i_spin) = & 
                                         ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
                                                 gradient_basis_basis(1:3),1)
  
     first_order_sigama(i_spin) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
                                        gradient_rho(1:3,i_spin,i_point),1)


!---------[1] non-explicit added ones for vibrations respons, same as E-field----- 
     first_order_H(i_basis,j_basis) =    &
     first_order_H(i_basis,j_basis) +    & 
        atomic_term * dir_tab_global(j_coord, j_atom, i_point)   +    & !(1) 
        point_term * v_hartree_gradient(i_point)                 +    & !(2)
        !----------(3.1) bassis_basis term-----------------------
        v2rho2(1,i_point)*basis_basis*first_order_rho(i_point)   + &
        v2rhosigma(1,i_point)*basis_basis*first_order_sigama(i_spin) + &
        !----------(3.2) gradient_basis_basis term---------------
        !----------(3.2.1) gradient_rho term--------------------
        2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)* & 
                                    first_order_rho(i_point) + &
        2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis(i_spin)*first_order_sigama(i_spin) + &
        !----------(3.2.2) first_order_gradient_rho term--------------------
        2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis(i_spin)

  end do ! n_points

  enddo
  enddo
 
!---------[2] explicit added ones for vibrations respons----- 
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1
 
     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)
 
       if(j_atom.eq.basis_atom(i_basis)) then 
       do i_point = 1, n_points, 1
         
         !---------first_order_gradient_basis_basis---------------- 
         !   d[ grad(X_u) X_v + X_u grad(X_v) ]/ d R_j
         ! = d[ grad(X_u) ]/dR_j X_v + grad(X_u) d(X_v)/dR_j + 
         !   d(X_u)/dR_j grad(X_v)   + X_u d[ grad(X_v) ]/dR_j
 
         do i_coord = 1,3 !loop over gradient coord
         first_order_gradient_basis_basis(i_coord) = partition_tab(i_point) * ( & 
                                 -hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point) & 
                                  *wave(j_compute,i_point)  &
                                 -gradient_basis_wave(i_compute,j_coord,i_point) & 
                                  *gradient_basis_wave(j_compute,i_coord,i_point)  )
         enddo
         gradient_rho_first_order_gradient_basis_basis(i_spin) = & 
                             ddot(3, gradient_rho(1:3,i_spin,i_point), 1, &
                                     first_order_gradient_basis_basis(1:3),1)
 
       first_order_H(i_basis,j_basis) =    &
       first_order_H(i_basis,j_basis) +    &
           partition_tab(i_point)* ( &  
           -gradient_basis_wave(i_compute,j_coord,i_point)*H_times_psi(j_compute,i_point)  & 
           -gradient_basis_wave(i_compute,j_coord,i_point)* &
           vrho(i_spin,i_point)*wave(j_compute,i_point) )   & 
            + 2.0d0*vsigma(1,i_point)*gradient_rho_first_order_gradient_basis_basis(i_spin) 
           
         if(flag_rel == REL_atomic_zora)then   
           first_order_H(i_basis,j_basis) =    &
           first_order_H(i_basis,j_basis) -    &
             partition_tab(i_point)* &
             wave(j_compute,i_point)*kinetic_gradient_basis_wave(i_compute,j_coord,i_point)
         endif
 
       enddo 
       endif 
 
       if(j_atom.eq.basis_atom(j_basis)) then 
       do i_point = 1, n_points, 1
 
         do i_coord = 1,3 !loop over gradient coord
         first_order_gradient_basis_basis(i_coord) = partition_tab(i_point) * ( & 
                                 -hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point) & 
                                  *wave(i_compute,i_point)  &
                                 -gradient_basis_wave(i_compute,i_coord,i_point) & 
                                  *gradient_basis_wave(j_compute,j_coord,i_point)  )
         enddo
         gradient_rho_first_order_gradient_basis_basis(i_spin) = & 
                             ddot(3, gradient_rho(1:3,i_spin,i_point), 1, &
                                     first_order_gradient_basis_basis(1:3),1)
 
       first_order_H(i_basis,j_basis) =    &
       first_order_H(i_basis,j_basis) +    &
           partition_tab(i_point)* ( &
           -gradient_basis_wave(j_compute,j_coord,i_point)*H_times_psi(i_compute,i_point)   &  
           -gradient_basis_wave(j_compute,j_coord,i_point)* & 
           vrho(i_spin,i_point)*wave(i_compute,i_point) )  & 
           + 2.0d0*vsigma(1,i_point)*gradient_rho_first_order_gradient_basis_basis(i_spin) 

         if(flag_rel == REL_atomic_zora)then   
           first_order_H(i_basis,j_basis) =    &
           first_order_H(i_basis,j_basis) -    &
             partition_tab(i_point)*&
             wave(i_compute,i_point)*kinetic_gradient_basis_wave(j_compute,j_coord,i_point)
         endif
 
       enddo 
       endif
  enddo
  enddo

  else 

  
    write(use_unit,'(2X,A)') "DFPT_reduce_memory is only support LDA, GGA, at present." 
    stop
  end if ! LDA, GGA 
end subroutine evaluate_first_order_H_reduce_memory
!******
