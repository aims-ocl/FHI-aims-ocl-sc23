!****s* FHI-aims/evaluate_first_order_rho_polarizability
!  NAME
!   evaluate_first_order_rho_polarizability
!  SYNOPSIS

subroutine evaluate_first_order_rho_polarizability( & 
           n_points,  & 
           KS_ev_compute,n_compute_c, i_basis_index, & 
           wave, &
           first_order_density_matrix, first_order_rho)

!  PURPOSE
!  Evaluates first_order_rho= rho_gradien+ first_order_DM*wave*wave
!  shanghui  2012.05.02

!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  use mpi_tasks, only: aims_stop, myid
  implicit none

!  ARGUMENTS

  integer :: n_points
  
  real*8, dimension(n_basis,n_states) :: KS_ev_compute

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave


  real*8:: first_order_density_matrix(3, n_basis, n_basis),Pone(n_compute_c,n_compute_c), tmp_rho(n_compute_c,n_points)
  real*8, dimension( 3, n_points) :: first_order_rho


! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix
! o wave -- basis functions
! o gradient_basis_wave -- gradient of basis functions
! o KS_ev_compute -- KS eigenvector

!
!  OUTPUT
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



 !     external functions
  real*8, external ::  ddot

 !  begin work
  first_order_rho(1:3,1:n_points) = 0.0d0

  do i_coord = 1, 3, 1 

    Pone = 0.d0
    tmp_rho = 0.d0

    do i_compute=1,n_compute_c
      i_basis=i_basis_index(i_compute)
      do j_compute=1,n_compute_c
        j_basis=i_basis_index(j_compute)         
        ! Replace 1st-order DM with n_compute_c dimensions for proper use with other arrays in dgemm
        Pone(i_compute,j_compute) = first_order_density_matrix(i_coord, i_basis, j_basis)  
      enddo
    enddo

    call dgemm("N","N",n_compute_c,n_points,n_compute_c,&
               1.d0,Pone,n_compute_c,&
               wave,n_max_compute_ham,0.d0,tmp_rho,n_compute_c)

     ! The following amounts to considering only the diagonal elements of the result of a matrix multiplication
     do i_point=1,n_points
       do i_compute=1,n_compute_c
         first_order_rho(i_coord,i_point)=first_order_rho(i_coord,i_point)+tmp_rho(i_compute,i_point)*wave(i_compute,i_point)
         !try with ddot instead?
       enddo
     enddo

  enddo


! old method: nested loops (slow)

!  first_order_rho(1:3,1:n_points) = 0.0d0
!
! !-------------second: first_order_rho---------------------------------
!  do i_coord = 1, 3, 1 
!
!       do i_point=1, n_points,1
!         
!          !-------shanghui add for safe---------------
!          first_order_rho(i_coord,  i_point)=0.0d0
!
!         do i_compute=1,n_compute_c
!         do j_compute=1,n_compute_c
!         
!          i_basis=i_basis_index(i_compute)
!          j_basis=i_basis_index(j_compute)          
!
!          first_order_rho(i_coord,  i_point)= &
!          first_order_rho(i_coord,  i_point)+ & 
!          first_order_density_matrix(i_coord, i_basis, j_basis) * & 
!          wave(i_compute,i_point)*wave(j_compute,i_point)           
!         enddo
!         enddo
!
!
!
!      end do  ! i_point
!  end do   ! i_coord


end subroutine evaluate_first_order_rho_polarizability
!---------------------------------------------------------------------
!******	 
