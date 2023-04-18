!****s* FHI-aims/evaluate_first_order_rho_moving_grid
!  NAME
!   evaluate_first_order_rho_moving_grid
!  SYNOPSIS

subroutine evaluate_first_order_rho_moving_grid( & 
           n_points,max_occ_number,  & 
           KS_ev_compute,n_compute_c, i_basis_index, & 
           wave,gradient_basis_wave, &
           first_order_density_matrix, first_order_rho_moving_grid)

!  PURPOSE
!  shanghui  2013.06.03

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
  real*8, dimension( 3, n_atoms, n_points) :: first_order_rho_moving_grid


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
  integer :: i_atom, j_atom

  real*8, dimension(3,n_atoms,max_occ_number, n_points) :: KS_orbital_gradient
  real*8, dimension(max_occ_number, n_points) :: KS_orbital

  real*8, dimension(3, n_atoms, n_points) :: rho_gradient_atom
  real*8, dimension( 3, n_atoms, n_points) :: first_order_DM_rho

 !     external functions
  real*8, external ::  ddot

  !  begin work
  rho_gradient_atom(1:3,1:n_atoms,1:n_points) = 0.0d0
  first_order_rho_moving_grid(1:3,1:n_atoms,1:n_points) = 0.0d0
  

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
                                 !shanghui change here for moving_grid--|
     KS_orbital_gradient(i_coord,basis_atom(i_basis),i_state,1:n_points)+ & 
     gradient_basis_wave(i_compute,i_coord,1:n_points)*                   &
     KS_ev_compute(i_basis,i_state)

    enddo
    enddo

  enddo

 do i_atom=1, n_atoms
     do i_point = 1, n_points, 1
        rho_gradient_atom(1:3,i_atom,i_point) =0.0d0
        do i_state=1, max_occ_number
           do j_atom=1, n_atoms

           if(j_atom.ne.i_atom) then
           rho_gradient_atom(1:3,i_atom,i_point) = &
           rho_gradient_atom(1:3,i_atom,i_point) + &
           2.0d0*KS_orbital(i_state,i_point)*       & ! here 2 is occ_nunber
           KS_orbital_gradient(1:3,j_atom,i_state,i_point)
           endif

           enddo

        enddo
     end do
 enddo

 
  rho_gradient_atom = 2.0d0 * rho_gradient_atom


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

      first_order_rho_moving_grid(i_coord, i_atom, i_point) =                     &
      rho_gradient_atom(i_coord, i_atom, i_point) +                    &
      first_order_DM_rho(i_coord, i_atom, i_point)


      end do  ! i_point
    end do ! i_atom
  end do   ! i_coord


end subroutine evaluate_first_order_rho_moving_grid
!---------------------------------------------------------------------
!******	 
