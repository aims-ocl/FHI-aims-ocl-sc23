!****s* FHI-aims/evaluate_first_order_rho
!  NAME
!   evaluate_first_order_rho
!  SYNOPSIS

subroutine evaluate_first_order_rho( & 
           n_points,max_occ_number,  & 
           KS_ev_compute,n_compute_c, i_basis_index, & 
           wave,gradient_basis_wave, &
           first_order_density_matrix, first_order_DM_rho, first_order_rho)

!  PURPOSE
!  Evaluates first_order_rho= rho_gradien+ first_order_DM*wave*wave
!  shanghui  2012.05.02

!  Now I add output also first_order_DM_rho
!  shanghui 2012.05.28
!  debug time we can refer to evaluate_density_gradient_denmat.f90
!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: max_occ_number
  
  real*8, dimension(n_basis,n_states) :: KS_ev_compute

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave


  real*8:: first_order_density_matrix(3, n_atoms, n_basis, n_basis)
  real*8, dimension( 3, n_atoms, n_points) :: first_order_DM_rho
  real*8, dimension( 3, n_atoms, n_points) :: first_order_rho


! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix
! o wave -- basis functions
! o gradient_basis_wave -- gradient of basis functions
! o KS_ev_compute -- KS eigenvector

!
!  OUTPUT
! o first_order_DM_rho : first_order_DM*wave*wave
! o first_order_rho : first_order_rho= rho_gradien+ first_order_DM*wave*wave
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



  !     counters
  integer :: i_coord
  integer :: i_point
  integer :: i_state 
  integer :: i_basis,j_basis,i_compute,j_compute
  integer :: i_atom

  real*8, dimension(3,n_atoms,max_occ_number, n_points) :: KS_orbital_gradient
  real*8, dimension(max_occ_number, n_points) :: KS_orbital

  real*8, dimension(3, n_atoms, n_points) :: rho_gradient_atom
!  real*8:: work(n_basis, n_points)
!  real*8:: temp_first_order(n_basis,n_basis) 

 !     external functions
  real*8, external ::  ddot

  !  begin work
  rho_gradient_atom(1:3,1:n_atoms,1:n_points) = 0.0d0
  first_order_DM_rho(1:3,1:n_atoms,1:n_points) = 0.0d0
  first_order_rho(1:3,1:n_atoms,1:n_points) = 0.0d0
  

  !--------------------first: rho_gradient_atom----------------------
  KS_orbital(1:max_occ_number,1:n_points)=0.0d0
  KS_orbital_gradient(1:3,1:n_atoms,1:max_occ_number,1:n_points)=0.0d0


 do i_state=1, max_occ_number 
     do i_compute=1,n_compute_c

     i_basis=i_basis_index(i_compute)

     KS_orbital(i_state,1:n_points)= &
     KS_orbital(i_state,1:n_points)+ &
     wave(i_compute,1:n_points)*       &
     KS_ev_compute(i_basis,i_state) 

     enddo
 enddo


 do i_state=1, max_occ_number 
  do i_coord = 1, 3, 1 
     do i_compute=1,n_compute_c

     i_basis=i_basis_index(i_compute)

     KS_orbital_gradient(i_coord,basis_atom(i_basis),i_state,1:n_points)= &
     KS_orbital_gradient(i_coord,basis_atom(i_basis),i_state,1:n_points)- &
     gradient_basis_wave(i_compute,i_coord,1:n_points)*                   &
     KS_ev_compute(i_basis,i_state)

     enddo
   enddo
  enddo

 do i_atom=1, n_atoms
     do i_point = 1, n_points, 1
        rho_gradient_atom(1:3,i_atom,i_point) =0.0d0
       do i_state=1, max_occ_number
        rho_gradient_atom(1:3,i_atom,i_point) = &
             rho_gradient_atom(1:3,i_atom,i_point) + &
             2.0d0*KS_orbital(i_state,i_point)*       & ! here 2 is occ_nunber
             KS_orbital_gradient(1:3,i_atom,i_state,i_point)
         enddo
     end do
 enddo



 
  rho_gradient_atom = 2.0d0 * rho_gradient_atom

!!!-------------test rho_gradient_atom
!!!    do i_basis=1,n_basis 
!!!    temp_first_order(i_basis,1:n_basis)= 2.0d0*  & 
!!!    KS_ev_compute(1,i_basis)*KS_ev_compute(1,i_basis)
!!!    enddo
!!!   
!!!      
!!!    CALL dgemm("N","N", n_basis, n_points, n_basis, &
!!!                    1.0d0, temp_first_order, n_basis,wave, & 
!!!                    n_basis,&
!!!                    0.0d0, work, n_basis )
!!!   
!!!    do i_point=1,n_points
!!!      
!!!    rho_gradient_atom(1,1,i_point)=   &
!!!    -work(1,i_point)*gradient_basis_wave(1,1, i_point)*2.0d0
!!!    enddo
!!!   

!!!  do i_point=1,n_points
!!!     write(use_unit,*) 'in evaluate_first_order_rho',i_point,rho_gradient_atom(1,1,i_point)
!!!  enddo
!!!  stop


  !-------------second: first_order_rho---------------------------------
  do i_coord = 1, 3, 1 
    do i_atom = 1, n_atoms, 1

       do i_point=1, n_points,1
         
          !-------shanghui add for safe---------------
          first_order_DM_rho(i_coord, i_atom, i_point)=0.0d0

         do i_compute=1,n_compute_c
         do j_compute=1,n_compute_c
         
          i_basis=i_basis_index(i_compute)
          j_basis=i_basis_index(j_compute)          

          first_order_DM_rho(i_coord, i_atom, i_point)= &
          first_order_DM_rho(i_coord, i_atom, i_point)+ & 
          first_order_density_matrix(i_coord, i_atom, i_basis, j_basis) * & 
          wave(i_compute,i_point)*wave(j_compute,i_point)           
         enddo
         enddo

      first_order_rho(i_coord, i_atom, i_point) =                     &
      rho_gradient_atom(i_coord, i_atom, i_point) +                    &
      first_order_DM_rho(i_coord, i_atom, i_point)

      !-------shanghui test v_gradient------------ 
      !first_order_DM_rho(i_coord, i_atom, i_point) =   & 
      !rho_gradient_atom(i_coord, i_atom, i_point) 
      !-------shanghui end test v_gradient---------

      end do  ! i_point
    end do ! i_atom
  end do   ! i_coord


end subroutine evaluate_first_order_rho
!---------------------------------------------------------------------
!******	 
